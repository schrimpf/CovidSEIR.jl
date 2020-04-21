"""
Functions to estimate a multiple region SEIR model.
"""
module MultiRegion

using ModelingToolkit, DifferentialEquations, Turing
#using CovidSEIR
using LinearAlgebra: dot
import ..RegionsData

"""
    Transfrom ODE variables to/from vector
"""
systemvec(S, E, I, C, CC, R, Rc, X) = [S, E, I, C..., CC, R, Rc, X]

"""
    Transfrom ODE variables to/from vector
"""
systemvars(vec) = (S=vec[1], E=vec[2], I=vec[3], C=vec[4:5],CC=vec[6], R=vec[7], Rc=vec[8], X=vec[9])

"""
    Transfrom ODE parameters to/from vector
"""
paramvec(β, γ, p, τ, a) = [β..., γ..., p..., τ, a]
paramvec(pm::NamedTuple) = paramvec(pm.β, pm.γ, pm.p, pm.τ, pm.a)

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
  return((β=β, γ=γ, p=p, τ=τ, a=a))
end

"""
   odeSEIR()

Sets up ODE for SEIR model with unconfirmed cases.

Returns an ODEProblem
"""
function odeSEIR()
  @parameters t β[1:3](t) γ[1:2] p[1:2] τ a
  @variables S(t) E(t) I(t) C[1:2](t) CC(t) R(t) Rc(t) X(t)
  @derivatives D'~t

  eqs = [D(S) ~ -S*( I*β[1]+ C[1]*β[2] + C[2]*β[3] ),
         D(E) ~ -a*E + S*( I*β[1] + C[1]*β[2] + C[2]*β[3] ),
         D(I) ~ a*E - γ[1]*I - p[1]*I - τ*I,
         D(C[1]) ~ τ*I - γ[1]*C[1] - p[1]*C[1],
         D(C[2]) ~ p[1]*(I + C[1]) - γ[2]*C[2] - p[2]*C[2],
         D(CC) ~ τ*I + p[1]*I,
         D(R) ~ γ[1]*I,
         D(Rc) ~ dot(γ,C),
         D(X) ~ p[2]*(C[2])
         ]

  de = ODESystem(eqs)
  f = ODEFunction(de, systemvec(S,E,I,C,CC,R,Rc,X), paramvec(β, γ, p, τ, a))
  prob = ODEProblem(f, systemvec(100.,1.,0.,zeros(2),0.,0.,0.,0.), (0.0, 100.0),
                    paramvec(fill(t->0, 3), zeros(2), zeros(2), 0., 0.))
  return(prob)
end

defaultpriors(;nx=2) = Dict(
  "a" => truncated(Normal(1/5, 3), 1/14, 1.0), # 1/incubation period
  "p[1]" => truncated(Normal(0.05, 0.3), 0, 1),
  "p[2]" => truncated(Normal(0.05, 0.3), 0, 1),
  "γ[1]" => truncated(Normal(0.133, 0.5), 0, 3),
  "γ[2]" => truncated(Normal(0.05, 0.3), 0, 1),
  "β[1]" => truncated(Normal(0.5, 1), 0, 10),
  "β[2]" => truncated(Normal(0.5, 1), 0, 10),
  "β[3]" => truncated(Normal(0.5, 1), 0, 10),
  "α" => MvNormal(zeros(nx), 1),
  "τ" => truncated(Normal(0.2, 2), 0, 10),
  "pE0" => truncated(Normal(0.001, 0.1), 0, 1),
  "sigD" => InverseGamma(2,3),
  "sigC" => InverseGamma(2,3),
  "sigRc" => InverseGamma(2,3),
  "sigS" => InverseGamma(2,3)
)


"""
     regionsmodel(data::RegionsData, priors=defaultpriors(),
                  ::Type{R}=Float64) where {R <: Real} = begni

Returns Turing model for multiple regions
"""
function regionsmodel(data::RegionsData, priors=defaultpriors(nx=size(data.X,1)))
  ode = odeSEIR()
  return(turingmodel(data.population, data.time, data.dead, data.active, data.recovered,
                     data.severe, data.X, ode, priors))
end

@model turingmodel(population, time, dead, active, recovered, severe, X,
                   ode, priors, ::Type{R}=Float64) where {R <: Real} =
begin

  # Allocate if simulating
  if dead===missing
      dead = Matrix{R}(undef, size(time))
  end
  if recovered===missing
    recovered = Matrix{R}(undef, size(time))
  end
  if active===missing
      active = Matrix{R}(undef, size(time))
  end
  if active===missing
    active = Matrix{R}(undef, size(time))
  end
  if severe===missing
    severe= Matrix{R}(undef, size(time))
  end

  # Priors
  L  = 2
  β = Vector{R}(undef, L+1)
  p = Vector{R}(undef, L)
  γ = Vector{R}(undef, L)
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
  end
  for i in 1:3
    β[i] ~ priors["β[$i]"]
  end
  α ~ priors["α"]
  τ ~ priors["τ"]
  pE0 ~ priors["pE0"]
  sigD  ~ priors["sigD"]
  sigC ~ priors["sigC"]
  sigRc ~ priors["sigRc"]
  sigS ~ priors["sigS"]
  N = length(population)
  eβ ~ MvNormal(zeros(N), 1)
  eτ ~ MvNormal(zeros(N), 1)
  ep ~ MvNormal(zeros(N), 1)

  # setup ODE
  for i in eachindex(population)
    B = β./population[i]
    βfn = [t->(B[1]*exp(α*X[:,Int64(round(t)),i] + eβ[i])),
           t->(B[2]*exp(α*X[:,Int64(round(t)),i] + eβ[i])),
           t->(B[3]*exp(eβ[i])) ]
    τi = τ*exp(eτ[i])
    pE0i = pE0*exp(ep[i])
    E0 = pE0i*population[i]
    S0 =population[i]*(1-pE0i)
    I0 = 0.0
    C0 = zeros(2)
    CC0= 0.0
    R0 = 0.0
    X0 = 0.0
    Rc0 = 0.0
    u0 = systemvec(S0,E0, I0, C0, CC0,R0, Rc0,X0)
    param = paramvec(βfn, γ, p, τi, a, ρ)
    probc = remake(ode, tspan = (0.0, Float64(maximum(time[:,i]))),u0=u0, p=param)
    sol = concrete_solve(probc, Tsit5(), u0, param; saveat=time[:,i])

    lu = size(sol,1)
    for t in 1:length(time[:,i])
      dead[t,i] ~ Normal(sol[lu,t], sigD)
      recovered[t,i] ~ Normal(sol[lu-1,t], sigRc)
      cumcases[t,i] ~ Normal(sol[6,t], sigC)
      severe[t,i] ~ Normal(sol[5,t], sigS)
    end

    return(active=active, dead=dead, recovered=recovered, severe=severe, E0=E0, sol=sol,
           β=β, γ=γ, τ=τ, p=p, a=a, sigD=sigD, sigRc=sigRc, sigC=sigC, sigS=sigS)
  end

end



end # Module
