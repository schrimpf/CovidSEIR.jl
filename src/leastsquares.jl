module LeastSquares

using DifferentialEquations, TransformVariables, Plots
import ColorSchemes
#using CovidSEIR
using LinearAlgebra: dot
import LineSearches
import Optim
import ..CountryData
import ..RegionsData
using Query, VegaLite, DataFrames


function mse(ode, data::CountryData)
  obj = let ode=ode, data=data
    function(θ::Vector{T}) where T
      if (length(θ) != length(ode.p) + 1)
        error("θ must be length $(length(ode.p) + 1)")
      end
      param = θ[1:(end-1)]
      pE0 = θ[end]
      u0 = [(1-pE0)*data.population, pE0*data.population, zeros(length(ode.u0)-2)...]
      prob = remake(ode, tspan=(0.0, maximum(data.time)), u0=u0, p=param)
      sol = solve(prob)
      se = zero(eltype(θ))
      lu = length(u0)
      for i in 1:length(data.time)
        se += (data.dead[i] - sol(data.time[i])[lu])^2 +
          (data.recovered[i] - sol(data.time[i])[lu-1])^2 +
          (data.active[i] - sum(sol(data.time[i])[4:5]))^2
      end
      se /= length(data.time)
    end
  end
end

function mse(ode, data::RegionsData, unpack)
  obj = let ode=ode, data=data
    function(θ::Vector{T}) where T
      (β, p, γ, τ, a, pE0, α, eβ, eτ, ep) = unpack(θ)
      N = length(data.population)
      uβ = Vector{eltype(eβ)}(undef, N)
      uβ[1:(N-1)] .= eβ
      uβ[N] = -sum(eβ)
      uτ = Vector{eltype(eτ)}(undef, N)
      uτ[1:(N-1)] .= eτ
      uτ[N] = -sum(eτ)
      up = Vector{eltype(ep)}(undef, N)
      up[1:(N-1)] .= ep
      up[N] = -sum(ep)
      se = zero(eltype(θ))
      for i in eachindex(data.population)
        B = β./data.population[i]
        βfn = [t->(B[1]*exp(dot(α,data.X[:,max(1,Int64(round(t))),i]) + uβ[i])),
               t->(B[2]*exp(dot(α,data.X[:,max(1,Int64(round(t))),i]) + uβ[i])),
               t->(B[3]*exp(uβ[i])) ]
        τi = τ*exp(uτ[i])
        pE0i = pE0*exp(up[i])
        E0 = pE0i*data.population[i]
        S0 =data.population[i]*(1-pE0i)
        I0 = 0.0
        C0 = zeros(2)
        CC0= 0.0
        R0 = 0.0
        X0 = 0.0
        Rc0 = 0.0
        u0 = [S0,E0, I0, C0..., CC0,R0, Rc0,X0]
        param = [βfn..., γ..., p..., τi, a]
        probc = remake(ode, tspan = (0.0, Float64(maximum(data.time[:,i]))),u0=u0, p=param)
        sol = concrete_solve(probc, Tsit5(), u0, param; saveat=data.time[:,i])
        if size(sol,2) < size(data.time,1)
          return(typeof(se)(Inf))
        end

        #sol = solve(probc)
        lu = size(sol,1)
        for t in 1:length(data.time[:,i])
          ismissing(data.dead[t,i])      || (se += (data.dead[t,i] - sol[lu,t])^2)
          ismissing(data.recovered[t,i]) || (se += (data.recovered[t,i] - sol[lu-1,t])^2)
          ismissing(data.cumcases[t,i])  || (se += (data.cumcases[t,i] - sol[6,t])^2)
          ismissing(data.severe[t,i])    || (se += (data.severe[t,i] - sol[5,t])^2)
        end
      end
      se /= length(data.time)
      return(se)
    end
  end
end

function leastsquares(data::RegionsData, ode)
  lbfgs = Optim.LBFGS(alphaguess=LineSearches.InitialStatic(alpha=1, scaled=true),
                      linesearch=LineSearches.BackTracking())
  trans = as( (β=as(Array, asℝ, 3),
               γ=as(Array, asℝ, 2),
               p=as(Array, asℝ, 2),
               τ=asℝ, a=asℝ,  pE0=asℝ,
               α=as(Array, asℝ, size(data.X,1)),
               eβ=as(Array, asℝ, length(data.population)-1),
               eτ=as(Array, asℝ, length(data.population)-1),
               ep=as(Array, asℝ, length(data.population)-1)) )
  obj = mse(ode, data, θ->transform(trans,θ))
  N = length(data.population)-1
  init = (β = [1.010040461970725e-6, 1.0000043030898965e-6, 0.7219713706236001], γ = [0.14124265615041995, 0.00857716350845789], p = [1.0012584545622381e-6, 0.04182070947753475], τ = 2.166207038767631, a = 0.035758335450852594, pE0 = 9.04616346465834e-6, α = [-0.101315343823542, -3.1263151170220085, -0.8364618759182466, -2.3339714749861376, -3.6307065120622872, -0.0014738347172689188, -0.8981184613104812, -4.5255211556090415e-5, -4.992148692652522],
          #β = [0.2, 0.1, 0.05],
          #γ = [0.05, 0.05],
          #p = [0.05, 0.05],
          #τ = 0.2, a = 1/5, pE0 = 1e-6,
          #α = -0.1*ones(size(data.X,1)),
          eβ = zeros(N), eτ=zeros(N), ep=zeros(N))
  x0 = inverse(trans, init)
  lb = (β = ones(3)*1e-8,
        γ = ones(2)*1e-8,
        p = ones(2)*1e-8,
        τ = 1e-6, a=1e-8, pE0=1/maximum(data.population),
        α = -5*ones(size(data.X,1)),
        eβ = -1*ones(N), eτ=-1*ones(N), ep=-1*ones(N) )
  ub = (β = ones(3)*3,
        γ = ones(2)*3,
        p = ones(2)*3,
        τ = 3, a=3, pE0=0.5,
        α = 0*ones(size(data.X,1)),
        eβ = 5*ones(N), eτ=5*ones(N), ep=5*ones(N) )
  nfree = length(x0) #-3*N
  opt = Optim.optimize(obj,
                       #x->obj(vcat(x, zeros(length(x0)-length(x)))),
                       inverse(trans,lb)[1:nfree],
                       inverse(trans,ub)[1:nfree],
                       x0[1:nfree],
                       Optim.Fminbox(lbfgs), Optim.Options(outer_iterations = 100,
                                                           show_trace=true,
                                                           extended_trace=false,
                                                           iterations=100),
                       autodiff=:forward)
  return(params=transform(trans,vcat(opt.minimizer, zeros(length(x0)-nfree))), opt=opt, obj=obj, x0=x0, trans=trans)
  #return(obj=obj, trans=trans, lb=lb, ub=ub, x0=x0, init=init)
end

function plotfits(data::RegionsData, ode, θ, trans)
  (β, p, γ, τ, a, pE0, α, eβ, eτ, ep) = transform(trans, θ)
  colors=ColorSchemes.colorschemes[:Set1_9]
  N = length(data.population)
  uβ = Vector{eltype(eβ)}(undef, N)
  uβ[1:(N-1)] .= eβ
  uβ[N] = -sum(eβ)
  uτ = Vector{eltype(eτ)}(undef, N)
  uτ[1:(N-1)] .= eτ
  uτ[N] = -sum(eτ)
  up = Vector{eltype(ep)}(undef, N)
  up[1:(N-1)] .= ep
  up[N] = -sum(ep)
  figs=[]
  sols=[]
  tsave = collect(1:(maximum(data.time)+7))
  Ts = maximum(tsave)
  for i in eachindex(data.population)
    B = β./data.population[i]
    (nx, T, N) = size(data.X)
    X = Array{eltype(data.X), 3}(undef, nx, Ts, N)
    X[:,1:T,:] .= data.X
    for t in (T+1):Ts
      X[:,t,:] .= X[:,t-1,:]
    end
    βfn = [t->(B[1]*exp(dot(α,X[:,max(1,Int64(round(t))),i]) + uβ[i])),
           t->(B[2]*exp(dot(α,X[:,max(1,Int64(round(t))),i]) + uβ[i])),
           t->(B[3]*exp(uβ[i])) ]
    τi = τ*exp(uτ[i])
    pE0i = pE0*exp(up[i])
    E0 = pE0i*data.population[i]
    S0 =data.population[i]*(1-pE0i)
    I0 = 0.0
    C0 = zeros(2)
    CC0= 0.0
    R0 = 0.0
    X0 = 0.0
    Rc0 = 0.0
    u0 = [S0,E0, I0, C0..., CC0,R0, Rc0,X0]
    param = [βfn..., γ..., p..., τi, a]
    probc = remake(ode, tspan = (0.0, Float64(maximum(tsave))),u0=u0, p=param)
    sol = concrete_solve(probc, Tsit5(), u0, param; saveat=tsave)
    lu = size(sol, 1)
    @show size(sol), size(tsave)
    r = data.recovered[:,i]
    if (all(ismissing.(r)))
      r[1] = 0
    end
    s = data.severe[:,i]
    if (all(ismissing.(s)))
      s[1] = 0
    end
    fig = plot(data.time[:,i], [data.dead[:,i], data.cumcases[:,i], s, r], color=colors[1:4]', label=["Dead" "cumulative cases" "severe"  "recovered"], title=data.name[i], legend=:topleft, linewidth=2)
    plot!(fig, tsave, [sol[lu,:], sol[6,:], sol[5,:], sol[lu-1,:]], color=colors[1:4]',
          label=nothing, linestyle=:dash, linewidth=1.5)
    push!(figs, fig)
    push!(sols, sol)
  end
  return(figs=figs, sol=sols, t=tsave)
end

function vegaplots(dat, f)

  labels = ["Dead","Cumulative Cases","Severe","Recovered"]
  lu = 9
  si = [lu, 6, 5, lu-1]
  di = [:dead, :cumcases, :severe, :recovered]
  foo = DataFrame(i=1:length(dat.name), pop=dat.population)
  foo = sort(foo, :pop, rev=true)

  P = 51
  obs = vcat([DataFrame(state=dat.name[i], group=".observed", var=lbl,
                        day = dat.time[:,i],
                        count=getfield(dat,d)[:,i]) for (lbl, d) in zip(labels, di) for i in foo.i[1:P]]...)
  est = vcat([DataFrame(state=dat.name[i], group="estimated", var=lbl,
                        day = f.t,
                        count=f.sol[i][s,:]) for (lbl, s) in zip(labels, si) for i in foo.i[1:P]]...)
  fdf = vcat(est, obs)

  vps=[(fdf |> @filter(_.state==st) |> @vlplot(:line, x=:day, y=:count, color=:var, wrap=:state, strokeDash=:group) ) for st in unique(fdf.state)]

  return(reduce((x,y)->[x; y], vps))

end


end
