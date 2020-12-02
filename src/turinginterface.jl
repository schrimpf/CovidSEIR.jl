function createturingmodel(data::CountryData, ode, systemt, paramt,
                           priors::NamedTuple{P, T}) where {P, T}

end

macro tt(foo::NamedTuple{P,T} where {P,T}
  println(P)
  println(T)
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
