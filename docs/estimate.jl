fileprefix = match(r"([a-zA-Z]+)_.+",ARGS[1]).captures[1]
println(ARGS)
println(fileprefix)
using CovidSEIR, CovidData, Turing, JLD2, Dates, DataFrames
countrynames = Dict( "us" => "US",
                     "china" => "China",
                     "korea" => "Korea, South",
                     "canada" => "Canada",
                     "italy" => "Italy")
                     
jmddir = normpath(joinpath(dirname(Base.find_package("CovidSEIR")),"..","docs","jmd"))
covdf = covidjhudata()
shift =fileprefix == "china" ?  65 : 1
dayt0 =  Dates.Date("2020-01-22")  - Dates.Day(shift)
dat = CountryData(covdf, countrynames[fileprefix], shift)
mdl = CovidSEIR.TimeVarying.countrymodel(dat, CovidSEIR.TimeVarying.defaultpriors())
println("Sampling for $fileprefix , ie $(countrynames[fileprefix])")
cc = Turing.psample(mdl, NUTS(0.65), 5000, 4)
JLD2.@save "$(jmddir)/$(fileprefix)_tv_$(Dates.today()).jld2" cc dayt0
