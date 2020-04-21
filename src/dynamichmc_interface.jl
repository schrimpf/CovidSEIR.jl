module DHMC
import ..CountryData

using DynamicHMC, Random, LogDensityProblems, Parameters, TransformVariables, Distributions, ForwardDiff, DifferentialEquations

"""
   Struct to represent a CovidSEIR problem for use with dynamichmc
"""
struct SEIRdhmc{P, V, O, I}
  data::CountryData
  priors::NamedTuple{P,V}
  ode::O
  "index of data.dead, data.recovered, and data.cases in solve(ode)"
  odeidx::I
  odeparams::Function
end

function (m::SEIRdhmc)(param)
  # priors
  @unpack data, priors, ode, odeidx, odeparams = m
  logp = mapreduce(x-> isa(x[1],Distribution) ? logpdf(x[1],x[2])
                   : mapreduce(y->logpdf(y[1],y[2]),+, zip(x[1],x[2]))
                   , +, zip(priors, param))
  u0 = [(1-param.pE0)*data.population, param.pE0*data.population, zeros(length(ode.u0))...]
  odep = odeparams(param)
  probc = remake(ode, tspan=(0.0, maximum(data.time)), u0=u0, p=odep)
  sol = concrete_solve(probc, Tsit5(), u0, odep; saveat=data.time)
  if size(sol,2) < length(data.time)
    return(logp + typeof(logp)(-Inf))
  end
  logp += logpdf(MvNormal(sol[odeidx.dead, :], param.sigD), data.dead)
  logp += logpdf(MvNormal(sol[odeidx.recovered, :], param.sigRc), data.recovered)
  meanc = vec(sum(sol[odeidx.active,:], dims=1))
  logp += logpdf(MvNormal(meanc, param.sigC), data.active)
  return(logp)
end

function transform_from_priors(priors)
  expr = Expr(:tuple)
  for (i,k) in enumerate(keys(priors))
    if isa(priors[k], Distribution)
      if (length(priors[k])==1)
        s = support(priors[k])
        lb= isinf(s.lb) ? -∞ : s.lb
        ub= isinf(s.ub) ? ∞ : s.ub
        push!(expr.args, :($k = as(Real, $(lb), $(ub))))
      else
        dims = size(priors[k])
        push!(expr.args, :($k = as(Array, asℝ, $(dims)...)))
      end
    else
      dims = length(priors[k])
      s = support(priors[k][1])
      lb= isinf(s.lb) ? -∞ : s.lb
      ub= isinf(s.ub) ? ∞ : s.ub
      push!(expr.args, :($k = as(Array, as(Real, $(lb), $(ub)), $(dims)...)))
    end
  end
  t = as((eval(expr) ))
  return(t)
end

function dhmcproblem(data::CountryData, priors::NamedTuple, ode,
                     odeidx::NamedTuple, odeparam::Function;
                     transform=transform_from_priors(priors),
                     ADbackend=:ForwardDiff)
  model = SEIRdhmc(data, priors, ode, odeidx, odeparam)
  P = TransformedLogDensity(transform, model)
  ∇P = ADgradient(ADbackend, P)
  return(prob=∇P, trans=transform)
end


end
