using Documenter, CovidSEIR

makedocs(;
         modules=[CovidSEIR],
         format=Documenter.HTML(),
         pages=[
           "Home" => "index.md",
           "Covid SEIR Estimates" => "covid.md",
           "Function Reference" => "functions.md",
         ],
         repo="https://github.com/schrimpf/CovidSEIR.jl/blob/{commit}{path}#L{line}",
         sitename="CovidSEIR.jl",
         #doctest=false,
         authors="Paul Schrimpf <schrimpf@mail.ubc.ca>"
         #assets=String[],
)


deploydocs(repo="github.com/schrimpf/CovidSEIR.jl.git")
