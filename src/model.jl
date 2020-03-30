"""
    Transfrom ODE variables to/from vector
"""
systemvec(S, E, I, C, R, Rc, X) = [S, E, I, C..., R, Rc, X]

"""
    Transfrom ODE variables to/from vector
"""
systemvars(vec) = (S=vec[1], E=vec[2], I=vec[3], C=vec[4:5], R=vec[6], Rc=vec[7], X=vec[8])

"""
    Transfrom ODE parameters to/from vector
"""
paramvec(β, γ, p, τ, a) = [β..., γ..., p..., τ, a]

"""
    Transfrom ODE parameters to/from vector
"""
paramvars(vec) = (β=vec[1:2], γ=vec[3:4], p=vec[5:6], τ=vec[7], a=vec[8])

"""
   odeSEIR()

Sets up ODE for SEIR model with unconfirmed cases.

Returns an ODEProblem
"""
function odeSEIR()
  @parameters t β[1:2] γ[1:2] p[1:2] τ a
  @variables S(t) E(t) I(t) C[1:2](t) R(t) Rc(t) X(t)
  @derivatives D'~t

  eqs = [D(S) ~ -S*( β[1]*I + dot(β, C) ),
         D(E) ~ S*( β[1]*I + dot(β, C) ) - a*E,
         D(I) ~ a*E - γ[1]*I - p[1]*I - τ*I,
         D(C[1]) ~ τ*I - γ[1]*C[1] - p[1]*C[1],
         D(C[2]) ~ p[1]*(I + C[1]) - γ[2]*C[2] - p[2]*C[2],
         D(R) ~ γ[1]*I,
         D(Rc) ~ dot(γ,C),
         D(X) ~ p[2]*(C[2])]

  de = ODESystem(eqs)
  f = ODEFunction(de, systemvec(S,E,I,C,R,Rc,X), paramvec(β, γ, p, τ, a))
  prob = ODEProblem(f, systemvec(100.,1.,0.,zeros(2),0.,0.,0.), (0.0, 100.0),
                    paramvec(zeros(2), zeros(2), zeros(2), 0., 0.))
  return(prob)
end

struct CountryData{T <: Real, I<:Number}
  population::T
  time::Vector{I}
  dead::Vector{T}
  recovered::Vector{T}
  active::Vector{T}
end

function CountryData(df::DataFrames.AbstractDataFrame, country::String, tstart=1)
  ss = (df.Country .== country)
  dead, conf, reco = if sum(ss .& .!ismissing.(df.Province))>0
    inc = ss .& .!ismissing.(df.Province)
    dead=DataFrames.by(df[inc, :], [:Country, :Date], :deaths => sum).deaths_sum
    conf=DataFrames.by(df[inc, :], [:Country, :Date], :confirmed => sum).confirmed_sum
    reco=DataFrames.by(df[inc, :], [:Country, :Date], :recovered => sum).recovered_sum
    dead, conf, reco
  else
    fill(missing, sum(ss)), fill(missing, sum(ss)), fill(missing, sum(ss))
  end
  if any(ismissing.(dead))
    dead = df[ss .& ismissing.(df.Province),:deaths]
  end
  if any(ismissing.(conf))
    conf = df[ss .& ismissing.(df.Province),:confirmed]
  end
  if any(ismissing.(reco))
    reco = df[ss .& ismissing.(df.Province),:recovered]
  end
  @assert length(dead)==length(conf)
  @assert length(dead)==length(reco)
  tvals = collect(tstart:(tstart+length(dead)-1))
  pop = unique(df[ss,:cpop])
  @assert length(pop)==1
  return(CountryData(Float64(pop[1]), tvals, Float64.(dead),
                     Float64.(reco), Float64.(conf-dead-reco)))
end

defaultcountrypriors() = Dict(
  "a" => truncated(Normal(1/5, 3), 1/14, 1.0), # 1/incubation period
  "p[1]" => truncated(Normal(0.05, 0.3), 0, 1),
  "p[2]" => truncated(Normal(0.05, 0.3), 0, 1),
  "γ[1]" => truncated(Normal(0.133, 0.5), 0, 3),
  "γ[2]" => truncated(Normal(0.05, 0.3), 0, 1),
  "β[1]" => truncated(Normal(0.5, 1), 0, 10),
  "β[2]" => truncated(Normal(0.5, 1), 0, 10),
  "τ" => truncated(Normal(0.2, 2), 0, 10),
  "pE0" => truncated(Normal(0.001, 0.1), 0, 1),
  "sigD" => InverseGamma(2,3),
  "sigC" => InverseGamma(2,3),
  "sigRc" => InverseGamma(2,3))  

"""
     countrymodel(data::CountryData, priors=defaultcountrypriors(),
                  ::Type{R}=Float64) where {R <: Real} = begni
    
Returns Turing model for single country
"""
function countrymodel(data::CountryData, priors=defaultcountrypriors())
  ode = odeSEIR()
  return(turingmodel1(data.population, data.time, data.dead, data.active, data.recovered, ode, priors))
end

@model turingmodel1(population, time, dead, active, recovered, ode,
                    priors=defaultcountrypriors(), ::Type{R}=Float64) where {R <: Real} =
begin
  # Priors
  L  = 2
  β = Vector{R}(undef, L)
  p = Vector{R}(undef, L)
  γ = Vector{R}(undef, L)
  def = defaultcountrypriors()
  for k in keys(def)
    if !(k in keys(priors))
      @warn "Key $k not found in priors. Using default value of $(def[k])"
      priors[k] = def[k]
    end
  end
  a ~ priors["a"] 
  for i in 1:L
    p[i] ~ priors["p[$i]"]
    γ[i] ~ priors["γ[$i]"]
    β[i] ~ priors["β[$i]"]
  end
  τ ~ priors["τ"]
  pE0 ~ priors["pE0"]
  sigD  ~ priors["sigD"]
  sigC ~ priors["sigC"]
  sigRc ~ priors["sigRc"]

  # setup ODE
  E0 = pE0*population
  S0 =population*(1-pE0)
  I0 = 0.0
  C0 = zeros(2)
  R0 = 0.0
  X0 = 0.0
  Rc0 = 0.0
  #[S, E, I..., C..., R, X, Rc], [β..., γ..., p..., τ..., a])
  u0 = systemvec(S0,E0, I0, C0, R0, X0, Rc0)
  param = paramvec(β./population, γ, p, τ, a)
  probc = remake(ode, tspan = (0.0, Float64(maximum(time))),u0=u0, p=param)
  sol = concrete_solve(probc, Tsit5(), u0, param; saveat=time) 

  if dead===missing
    dead = Vector{R}(undef, length(time))
  end
  if recovered===missing
    recovered = Vector{R}(undef, length(time))
  end
  if active===missing
    active = Vector{R}(undef, length(time))
  end
  lu = length(u0)
  for i in 1:length(time)
    t = time[i]
    #S, E, I, C, R, Rc, X = systemvars(sol[:,i])
    dead[i] ~ Normal(sol[lu,i], sigD)
    recovered[i] ~ Normal(sol[lu-1,i], sigRc)
    active[i] ~ Normal(sum(sol[4:5,i]), sigC)
  end
  return(active=active, dead=dead, recovered=recovered,E0=E0, sol=sol,
         β=β, γ=γ, τ=τ, p=p, a=a, sigD=sigD, sigRc=sigRc, sigC=sigC)
end

"""
    priorreport((priors=defaultcountrypriors(), T=100, population=1e7)

Create tables and figures summarizing priors.
"""
function priorreport(priors=defaultcountrypriors(), T=100, population=1e7;
                     colors=ColorSchemes.colorschemes[:Set1_9])
  sim = turingmodel1(population, 1:T, missing, missing, missing, odeSEIR(), priors)
  d = sim()
  alpha = 0.05
  toplot(d) = [d.active, d.recovered, d.dead]./population
  all = toplot(d)
  fig = Vector{typeof(plot())}(undef, length(all))
  draws = 1000
  for i in 1:length(all)
    fig[i] = plot(1:T, toplot(d)[i], alpha=alpha, legend=false, color=colors[i])                
  end
  for s in 2:draws
  d = sim()
    all = hcat.(all, toplot(d))
    for i in 1:length(all)
      fig[i]=plot!(fig[i], 1:T, toplot(d)[i], alpha=alpha, legend=false, color=colors[i])
    end
  end
  
  q(p) = map(d->mapslices(x->quantile(x, p),d, dims=2), all)
  lo =q(0.05)
  hi = q(0.95)
  
  for i in 1:length(all)
    fig[i]=plot!(fig[i], 1:T, mean(all[i], dims=2), alpha=1, 
                 color=colors[i], linewidth=3, linecolor=:black, ribbon=(mean(all[i],dims=2).-lo[i],
                                                       hi[i] .- mean(all[i],dims=2)), fillalpha=0.3)
  end
  plot!(fig[3], xlabel="Days")
  plot!(fig[2], ylabel="Portion of population")
  
  combofig=plot(fig..., title=["Active Confirmed Cases" "Recoveries" "Deaths"], link=:x, layout=(3,1))

  @model pm(priors, ::Type{R}=Float64) where  {R <: Real} =begin
    L = 2
    β = Vector{R}(undef, 2)
    p = Vector{R}(undef, 2)
    γ = Vector{R}(undef, 2)
    a ~ priors["a"] 
    for i in 1:L
      p[i] ~ priors["p[$i]"]
      γ[i] ~ priors["γ[$i]"]
      β[i] ~ priors["β[$i]"]
    end
    τ ~ priors["τ"]
    pE0 ~ priors["pE0"]
    sigD  ~ priors["sigD"]
    sigC ~ priors["sigC"]
    sigRc ~ priors["sigRc"]
  end
  pc = sample(pm(defaultcountrypriors()), NUTS(0.65), 500)
  
  return(all=combofig, figs=fig, tbl=describe(pc))
end

"""
    simtrajectories(cc::AbstractMCMC.AbstractChains,
                         data::CountryData, ts;
                         ic=Iterators.product(StatsBase.sample(1:size(cc,1),300, replace=false),
                                              1:size(cc,3)))


Simulate trajectories based on parameter values in chain `cc`.
"""
function simtrajectories(cc::AbstractMCMC.AbstractChains,
                         data::CountryData, ts;
                         ic=Iterators.product(StatsBase.sample(1:size(cc,1),300, replace=false),
                                              1:size(cc,3)))
  ode = odeSEIR()
  df = DataFrames.DataFrame()
  for (i, c) in ic
    a = cc.value[i, "a", c]
    pE0 = cc.value[i, "pE0", c]
    p = [cc.value[i,"p[1]",c] cc.value[i,"p[2]",c]]
    β = [cc.value[i,"β[1]",c] cc.value[i,"β[2]",c]]
    γ = [cc.value[i,"γ[1]",c] cc.value[i,"γ[2]",c]]
    sigC = cc.value[i,"sigC",c]
    sigD = cc.value[i,"sigD",c]
    sigRc = cc.value[i,"sigRc",c]
    τ = cc.value[i,"τ",c]
    param=[(β./data.population)..., γ..., p..., τ, a]
    u0 = [data.population*(1-pE0), data.population*pE0, zeros(6)...]
    devars = [:S, :E, :I, :C1, :C2, :R, :Rc, :X]
    probc = remake(ode, tspan = (0.0, Float64(maximum(ts))),u0=u0, p=param)
    sol = solve(probc, Tsit5())
    X = hcat(sol(ts)...)'
    ndf = DataFrames.DataFrame(:iter=>i, :chain=>c, :t=>ts)
    for (v, x) in zip(devars, eachcol(X))
      ndf[!,v] = x
    end    
    n = size(ndf,1)
    ndf[!,:dead] = ndf[!,:X] .+ randn(n)*sigD
    ndf[!, :active] = ndf[!, :C1] .+ ndf[!, :C2] .+ randn(n)*sigC
    ndf[!,:recovered] = ndf[!,:Rc] .+ randn(n)*sigRc
    append!(df, ndf)
  end
  return(df)
end
function varmatrix(adf, vars, fn=x->x)
  vars = vec(vars)
  meanpred = Matrix(adf[:,Symbol.(vars.*"_mean")])
  q5= Matrix(adf[:,Symbol.(vars.*"_function")])
  q95= Matrix(adf[:,Symbol.(vars.*"_function_1")])
  return(fn.(meanpred), fn.(q5), fn.(q95))
end

"""
    plotvars(simdf::DataFrames.AbstractDataFrame,
             data::CountryData;
             dayt0=Dates.Date("2020-01-21"), # one day before JHU data begins
             colors=ColorSchemes.colorschemes[:Set1_9])

Create plots showing fit and implications of simulated trajectories.
"""
function plotvars(df::DataFrames.AbstractDataFrame,
                  data::CountryData;
                  dayt0=Dates.Date("2020-01-21"), # one day before JHU data begins
                  colors=ColorSchemes.colorschemes[:Set1_9])
  date = dayt0 .+ Dates.Day.(data.time)
  df[!,:cfr] = df.dead./(df.active+df.recovered+df.dead)
  df[!,:tcfr] = df.dead./(df.active+df.recovered+df.dead+df.Rc+df.I)
  df[!,:pundect] = df.I./(df.C1+df.C2+df.I)
  df[!,:pmild] = df.C1./(df.C1+df.C2)
  df[!,:tpmild] = (df.I+df.C1)./(df.C1+df.C2+df.I)  
  adf = DataFrames.aggregate(DataFrames.groupby(df, :t), [mean, std, x->quantile(x, 0.05), x->quantile(x, 0.95)])
  adf[!,:date] = dayt0 .+ Dates.Day.(adf.t)
  minval = 0.1 # so log scales don't get log(0)

  # fit of observables
  fit = plot(date, max.(minval, [data.dead data.recovered data.active]),
             label=["dead" "recovered" "active"],
             color=colors[1:3]', linestyle=:solid, linewidth=2)
  meanpred, q5, q95 = varmatrix(adf,["dead","recovered","active"], x->max.(x, minval))
  daylim = (Dates.today() - Dates.Day(21), Dates.today() + Dates.Day(14))
  fit=plot!(fit,adf.date, meanpred, ribbon=(meanpred-q5,q95-meanpred), fillalpha=0.2, label=nothing,
            linestyle=:dash, xlim=Dates.value.(daylim), color=colors[1:3]',
            xrotation=25,
            ylims=(minval, maximum(data.active)*1.2),           
            yscale=:identity)
  
  # plots of trajectories
  vars = [ :S => "susceptible",
           :E => "exposed",
           :I => "undetected mild cases",
           :C1 => "detected mild cases",
           :C2 => "severe cases",
           :active => "detected cases",
           :R => "undetected recoveries",
           :Rc => "detected recoveries",
           :X => "deaths",
           :cfr=> "detected CFR",
           :tcfr=>"true CFR",
           :pundect=>"P(infections undetected)",
           :pmild=>"P(mild|detected)",
           :tpmild=>"P(mild)" ]
  tfigs = Vector{typeof(fit)}(undef, 0)
  daylim = (Dates.today() - Dates.Day(21), Dates.today() + Dates.Day(60))
  c = 1
  for v in vars    
    meanpred, q5, q95 = varmatrix(adf,[String(v[1])], x->x)
    fig = plot(adf.date, meanpred, title=v[2], ribbon=(meanpred-q5, q95-meanpred),
               legend=false, xlim=Dates.value.(daylim), xrotation=25, color=colors[c])
    if v[1] in [:cfr, :tcfr]
      fig = plot!(fig, ylim=(0.0, 0.3))
    elseif v[1] in [:pundect, :pmild, :tpmild]
      fig = plot!(fig, ylim=(0.0, 1.0))
    end
    push!(tfigs, fig)
    c = mod(c, length(colors)) + 1
  end  
  return(fit=fit, trajectories=tfigs)
end
