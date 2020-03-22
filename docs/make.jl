using Documenter, CovidSEIR

makedocs(;
    modules=[CovidSEIR],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/schrimpf/CovidSEIR.jl/blob/{commit}{path}#L{line}",
    sitename="CovidSEIR.jl",
    authors="Paul Schrimpf <paul.schrimpf@gmail.com>",
    assets=String[],
)
