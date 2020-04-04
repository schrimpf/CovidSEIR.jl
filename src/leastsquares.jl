module LeastSquares

using DifferentialEquations
#using CovidSEIR
using LinearAlgebra: dot
import ..CountryData

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

end
