module TimeVarying

using ModelingToolkit, DifferentialEquations, Turing
#using CovidSEIR
using LinearAlgebra: dot
import ..CountryData

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
paramvec(β, γ, p, τ, a, ρ) = [β..., γ..., p..., τ, a, ρ...]

"""
    Transfrom ODE parameters to/from vector
"""
function paramvars(vec)
  β=vec[1:3]
  c = length(β)
  γ=vec[c .+ (1:2)]
  c += length(γ)
  p=vec[c .+ (1:2)]
  c += length(p)
  τ=vec[c+1]
  c += length(τ)
  a=vec[c+1]
  c += length(a)
  ρ=vec[c .+ (1:2)]
  return((β=β, γ=γ, p=p, τ=τ, a=a,ρ=ρ))
end

"""
   odeSEIR()

Sets up ODE for SEIR model with unconfirmed cases.

Returns an ODEProblem
"""
function odeSEIR()
  @parameters t β[1:3] γ[1:2] p[1:2] τ a ρ[1:2]
  @variables S(t) E(t) I(t) C[1:2](t) R(t) Rc(t) X(t)
  @derivatives D'~t
  ρ3 = -1

  eqs = [D(S) ~ -S*( (1-ρ[1]/(1+exp((t-ρ[2])*ρ3)))*(I*(β[1]+β[3]) + β[1]*C[1]) + C[2]*β[2] ),
         D(E) ~ -a*E + S*( (1-ρ[1]/(1+exp((t-ρ[2])*ρ3)))*(I*(β[1]+β[3]) + β[1]*C[1]) + C[2]*β[2] ),
         D(I) ~ a*E - γ[1]*I - p[1]*I - τ*I,
         D(C[1]) ~ τ*I - γ[1]*C[1] - p[1]*C[1],
         D(C[2]) ~ p[1]*(I + C[1]) - γ[2]*C[2] - p[2]*C[2],
         D(R) ~ γ[1]*I,
         D(Rc) ~ dot(γ,C),
         D(X) ~ p[2]*(C[2])]

  de = ODESystem(eqs)
  f = ODEFunction(de, systemvec(S,E,I,C,R,Rc,X), paramvec(β, γ, p, τ, a, ρ))
  prob = ODEProblem(f, systemvec(100.,1.,0.,zeros(2),0.,0.,0.), (0.0, 100.0),
                    paramvec(zeros(3), zeros(2), zeros(2), 0., 0.,zeros(2)))
  return(prob)
end

defaultpriors() = Dict(
  "a" => truncated(Normal(1/5, 3), 1/14, 1.0), # 1/incubation period
  "p[1]" => truncated(Normal(0.05, 0.3), 0, 1),
  "p[2]" => truncated(Normal(0.05, 0.3), 0, 1),
  "γ[1]" => truncated(Normal(0.133, 0.5), 0, 3),
  "γ[2]" => truncated(Normal(0.05, 0.3), 0, 1),
  "β[1]" => truncated(Normal(0.5, 1), 0, 10),
  "β[2]" => truncated(Normal(0.5, 1), 0, 10),
  "β[3]" => truncated(Normal(0.5, 1), 0, 10),
  "τ" => truncated(Normal(0.2, 2), 0, 10),
  "pE0" => truncated(Normal(0.001, 0.1), 0, 1),
  "sigD" => InverseGamma(2,3),
  "sigC" => InverseGamma(2,3),
  "sigRc" => InverseGamma(2,3),
  "ρ[1]" => truncated(Normal(0.5, 2), 0, 1), # reduction in β
  "ρ[2]" => truncated(Normal(30, 30), 0, 100) #, # time of reduction
  #"ρ[3]" => truncated(Normal(-0.5, 2), -Inf, 0) # speed of reduction
)
  

"""
     countrymodel(data::CountryData, priors=defaultcountrypriors(),
                  ::Type{R}=Float64) where {R <: Real} = begni
    
Returns Turing model for single country
"""
function countrymodel(data::CountryData, priors=defaultpriors())
  ode = odeSEIR()
  return(turingmodel1(data.population, data.time, data.dead, data.active, data.recovered, ode, priors))
end

@model turingmodel1(population, time, dead, active, recovered, ode,
                    priors=defaultpriors(), ::Type{R}=Float64) where {R <: Real} =
begin
  # Priors
  L  = 2
  β = Vector{R}(undef, L+1)
  p = Vector{R}(undef, L)
  γ = Vector{R}(undef, L)
  ρ = Vector{R}(undef, 2)
  def = defaultpriors()
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
    ρ[i] ~ priors["ρ[$i]"]
  end
  for i in 1:3
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
  u0 = systemvec(S0,E0, I0, C0, R0, Rc0,X0)
  param = paramvec(β./population, γ, p, τ, a, ρ)
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
    #t = time[i]
    #S, E, I, C, R, Rc, X = systemvars(sol[:,i])
    dead[i] ~ Normal(sol[lu,i], sigD)
    recovered[i] ~ Normal(sol[lu-1,i], sigRc)
    m = sum(sol[4:5,i])
    active[i] ~ Normal(m,sigC)
  end
  return(active=active, dead=dead, recovered=recovered,E0=E0, sol=sol,
         β=β, γ=γ, τ=τ, p=p, a=a, sigD=sigD, sigRc=sigRc, sigC=sigC)
end


end
