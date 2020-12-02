module CovidSEIR

import DataFrames
using CovidData
import Dates, AbstractMCMC, StatsBase,ColorSchemes
using DifferentialEquations: concrete_solve, solve, remake, Tsit5, ODEProblem
using ModelingToolkit: @parameters, @variables, @derivatives, ODESystem, ODEFunction, ~, Variable, Differential
using Distributions
using Turing
using LinearAlgebra: dot
using Statistics
using Plots, StatsPlots

include("data.jl")
include("model.jl")
include("timevarying.jl")
include("leastsquares.jl")
include("dynamichmc_interface.jl")
include("multiregion.jl")

export CountryData,
  countrymodel,
  priorreport,
  simtrajectories,
  plotvars

end # module
