module RT
using DynamicHMC, TransformVariables, LogDensityProblems, Parameters, Distributions
using LinearAlgebra, Plots
using StaticArrays

struct RtRW{C,P,V,T, T0}
  dlogk::C
  X::T
  X0::T0
  priors::NamedTuple{P,V}
  MA1::Bool
end

function RtRW(d::Array{Vector{T}}, p::NamedTuple{P,V}) where {T, P , V}
  RtRW(d, [ones(length(ds),1) for ds in d],
       fill(ones(1), length(d)), p, true)
end

RtRW(d, x, x0, p) = RtRW(d,x,x0,p,false)

function (m::RtRW)(param)
  @unpack dlogk, X, X0, priors = m
  @unpack σR, σk, σR0, γ, ρ, α, α0 = param
  logp = logpdf(priors.σR, σR) +
    logpdf(priors.γ, γ) +
    logpdf(priors.σR0, σR0) +
    #logpdf(priors.μR0, μR0) +
    logpdf(priors.σk, σk) +
    logpdf(priors.α, α) +
    logpdf(priors.α0, α0) +
    logpdf(priors.ρ,ρ)

  # Note that this not the usual Kalman filter because the noise in
  # dlogk is an MA(1) instead of independent. We use
  # https://arxiv.org/pdf/1909.10582.pdf
  # and adopt their notation
  A = m.MA1 ? SMatrix{2,2}([2*σk^2  -σk^2; -σk^2  2*σk^2]) : SMatrix{2,2}([2*σk^2 0; 0 2*σk^2])
  S = length(dlogk)
  T = maximum(length.(dlogk))
  meanR = Vector{typeof(σR)}(undef, T)
  zhat = Vector{typeof(σR)}(undef, T)
  varR, varRprior, K, varz, zcoef = ma1kalman_variance(T, σR0^2, ρ, σR^2, γ, A)
  for s in 1:S
    z = dlogk[s] .+ γ .- γ* (X[s]*α)
    μR0 = dot(X0[s],α0)
    ma1kalman_mean!(meanR, zhat, z, ρ, γ, μR0, K, zcoef) #, X[s]*α)
    for t in 1:length(dlogk[s])
      μ = zhat[t]
      σ = sqrt(varz[t])
      logp += logpdf(Normal(μ, σ), z[t])
    end
  end
  return(logp)
end

function ma1kalman_variance(T, V0, F, W, H, A)
  # Note that this not the usual Kalman filter because the noise in
  # dlogk is an MA(1) instead of independent. We use
  # https://arxiv.org/pdf/1909.10582.pdf
  # and adopt their notation
  varR = Vector{typeof(V0)}(undef,T)
  varRprior = Vector{typeof(V0)}(undef,T)
  K = Vector{typeof(V0)}(undef,T)
  varz = Vector{typeof(V0)}(undef, T)
  zcoef = Matrix{typeof(V0)}(undef,2, T)
  varZ = MMatrix{2,2,typeof(V0)}(zeros(4))
  covzZ = MVector{2, typeof(V0)}(zeros(2))
  for t in 1:T
    varZ .= 0
    if t == 1
      varRprior[t] = V0
    else
      varRprior[t] = F*varR[t-1]*F' + W
      varZ[2,2] = H*varRprior[t-1]*H'
    end
    varZ[1,1] = H*varRprior[t]*H'
    varZ .+= A
    varz[t] = varZ[1,1]
    covzZ[1] = varZ[1,2]
    ivarZ = inv(varZ)
    zcoef[:,t] .= vec(covzZ'*ivarZ)
    iL = 1/(varZ[1,1] + dot(zcoef[:,t], covzZ))
    J = varRprior[t]*H
    K[t] = J*iL
    varR[t] =  varRprior[t]-K[t]*J' # P_t
  end
  return(varR, varRprior, K, varz, zcoef)
end

function ma1kalman_mean!(x, zhat, z, F, H, x0, K, zcoef)
  T = length(z)
  elag = @MVector zeros(eltype(zcoef[1,1]*z[1]), 2)
  for t in 1:T
    xminus = F*(t==1 ? x0 : x[t-1])
    ez = H*xminus
    zhat[t] = ez
    zhat[t] += zcoef[:,t]'*elag
    x[t] = xminus + K[t]*(z[t] - zhat[t])
    elag[2] = elag[1]
    elag[1] = z[t] - H*xminus
  end
end

function ma1kalman_mean!(x, zhat, z, F, H, x0, K, zcoef, xshift)
  T = length(z)
  elag = @MVector zeros(eltype(zcoef[1,1]*z[1]), 2)
  for t in 1:T
    xminus = F*(t==1 ? x0 : x[t-1] ) + xshift[t]
    ez = H*xminus
    zhat[t] = ez
    #zhat[t] += zcoef[:,t]'*elag
    x[t] = xminus + K[t]*(z[t] - zhat[t])
    elag[2] = elag[1]
    elag[1] = z[t] - H*xminus
  end
end


function kalman(dlogk, σR, σk, σR0, μR0, γ, ρ=1, Xα=zeros(length(dlogk)); MA1=false)
  T = length(dlogk)
  A = MA1 ?  SMatrix{2,2}([2*σk^2  -σk^2; -σk^2  2*σk^2]) :  SMatrix{2,2}([2*σk^2  0; 0  2*σk^2])
  meanR = Vector{typeof(σR)}(undef, T)
  zhat = Vector{typeof(σR)}(undef, T)
  varR, varRprior, K, varz, zcoef = ma1kalman_variance(T, σR0^2, ρ, σR^2, γ, A)
  z = dlogk .+ γ .- γ* (Xα)
  ma1kalman_mean!(meanR, zhat, z, ρ, γ, μR0, K, zcoef)
  meanR .+= Xα
  zhat .+= -γ .+γ*(Xα)
  return(meanR, varR, zhat)
end


function smoother(x::AbstractVector; n=7, w=pdf(Normal(), range(-3, 3, length=n)))
  sx = Vector{Union{Missing, Float64}}(undef, length(x))
  n = length(w)
  if n % 2 != 1
    error("only odd length windows allowed")
  end
  shift = -(n÷2):(n÷2)
  for i in 1:length(sx)
    s = findfirst((i .+ shift) .> 0)
    l = findlast((i .+ shift) .< length(x))
    sx[i] = sum(w[s:l].*x[i .+ shift[s:l]])./sum(w[s:l])
  end
  return(sx)
end

function plotpostr(dates, dlogk, post, X, X0; Δ=1)
  k = [RT.kalman(dlogk, p.σR, p.σk, p.σR0, dot(X0,p.α0), p.γ, p.ρ, X*p.α) for p in post];
  γ = [p.γ for p in post];
  Xa = hcat([X*p.α for p in post]...)
  meanR = hcat([x[1] for x in k]...)./Δ;
  varR = hcat([x[2] for x in k]...)./Δ^2;
  zhat = hcat([x[3] for x in k]...);
  c = "black"
  figr = plot(dates, mean(meanR, dims=2), ribbon=1.64*mean(sqrt.(varR),dims=2), color=c, fillalpha=0.2,
              linewidth=1.5, legend=:none, ylab="Rₜ")
  r=([quantile(meanR[t,:] - 1.64*sqrt.(varR[t,:]), 0.05) for t in 1:size(meanR,1)],
     [quantile(meanR[t,:] + 1.64*sqrt.(varR[t,:]), 0.95) for t in 1:size(meanR,1)])
  figr = plot!(figr, dates, zeros(length(r[1])), ribbon=(-r[1], r[2]), color=c,
               linewidth=0, ylim=nothing, fillalpha=0.2)
  T = length(dlogk)
  rt = [mean(dlogk[t]./γ .+ 1) for t in 1:T]
  lo = [quantile(dlogk[t]./γ .+ 1, 0.05) for t in 1:T]
  hi = [quantile(dlogk[t]./γ .+ 1, 0.95) for t in 1:T]
  figr = scatter!(figr, dates, rt, yerror=(rt-lo, hi-rt))
  figr = hline!(figr, [1.], color="red", linewidth=1.5, linestyle=:dash)
  return(figr)
end



end
