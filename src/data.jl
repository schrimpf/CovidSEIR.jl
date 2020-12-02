################################################################################
"""
   struct CountryData{T <: Real, I<:Number}

Container for data for estimation.
"""
struct CountryData{T <: Real, I<:Number}
  population::T
  "time of observations"
  time::Vector{I}
  "dead[i] = number of deaths at time[i]"
  dead::Vector{T}
  "recovered[i] = number of deaths at time[i]"
  recovered::Vector{T}
  "active[i] = number of deaths at time[i]"
  active::Vector{T}
end

"""
    CountryData(df::DataFrames.AbstractDataFrame, country::String, tstart=1)

Construct `CountryData` from observation of `country` in `df`.

Will set the time variable in the ODE to be `tstart` at the first date in the data frame.
"""
function CountryData(df::DataFrames.AbstractDataFrame, country::String, tstart=1;
                     countryvar=:Country, datevar=:Date,
                     deadvar=:deaths, casevar=:confirmed,
                     recoveredvar=:recovered, popvar=:cpop)
  ss = (df.Country .== country)
  dead, conf, reco = if sum(ss .& .!ismissing.(df.Province))>0
    inc = ss .& .!ismissing.(df.Province)
    dead=DataFrames.by(df[inc, :], [countryvar, datevar], deadvar => sum).deaths_sum
    conf=DataFrames.by(df[inc, :], [countryvar, datevar], casevar => sum).confirmed_sum
    reco=DataFrames.by(df[inc, :], [countryvar, datevar], recoveredvar => sum).recovered_sum
    dead, conf, reco
  else
    fill(missing, sum(ss)), fill(missing, sum(ss)), fill(missing, sum(ss))
  end
  if any(ismissing.(dead))
    dead = df[ss .& ismissing.(df.Province),deadvar]
  end
  if any(ismissing.(conf))
    conf = df[ss .& ismissing.(df.Province),casevar]
  end
  if any(ismissing.(reco))
    reco = df[ss .& ismissing.(df.Province),recoveredvar]
  end
  @assert length(dead)==length(conf)
  @assert length(dead)==length(reco)
  tvals = collect(tstart:(tstart+length(dead)-1))
  pop = unique(df[ss,popvar])
  @assert length(pop)==1
  return(CountryData(Float64(pop[1]), tvals, Float64.(dead),
                     Float64.(reco), Float64.(conf-dead-reco)))
end


"""
   struct RegionsData{T <: Real, I<:Number}

Container for data for estimation.
"""
struct RegionsData{T <: Real, I<:Number}
  population::Vector{T}
  "time of observations"
  time::Matrix{I}
  "dead[t,i] = number of deaths at time[t,i]"
  dead::Matrix{Union{Missing,T}}
  "recovered[t,i] = number of recoveries at time[t,i]"
  recovered::Matrix{Union{Missing,T}}
  "cumcases[t,i] = number of cumulative cases at time[t,i]"
  cumcases::Matrix{Union{Missing,T}}
  "severe[t,i] = number of active severe cases at time[t,i]"
  severe::Matrix{Union{Missing,T}}
  "X[:,t,i] are covariates that shift infection rates"
  X::Array{T,3}
  name::Vector{String}
end

"""
    RegionsData(df::DataFrames.AbstractDataFrame, tstart=1)

Construct `CountryData` from observation of `country` in `df`.

Will set the time variable in the ODE to be `tstart` at the first date in the data frame.
"""
function RegionsData(df::DataFrames.AbstractDataFrame, tstart=1;
                     idvar=:fips, datevar=:date, deadvar=:deaths, casevar=:cases,
                     recoveredvar=Symbol("recovered.ctp"), severevar=:hospitalized,
                     popvar=Symbol("Population.2018"), namevar=:state,
                     policyvars=[Symbol("State.of.emergency"),
                                 Symbol("Date.closed.K.12.schools"),
                                 Symbol("Stay.at.home..shelter.in.place"),
                                 Symbol("Closed.restaurants.except.take.out"),
                                 Symbol("Closed.non.essential.businesses")],
                     otherxvars=[:workplaces_percent_change_from_baseline,
                                 :grocery_and_pharmacy_percent_change_from_baseline,
                                 :retail_and_recreation_percent_change_from_baseline,
                                 :percentchangebusinesses])


  ids = unique(df[!,idvar])
  dates = sort(unique(df[!,datevar]))
  N = length(ids)
  T = length(dates)
  nx = length(policyvars) + length(otherxvars)
  X=Array{Float64, 3}(undef, nx, T, N)
  dead = Matrix{Union{Missing, Float64}}(missing, T,N)
  tvals = Matrix{Int64}(undef, T, N)
  cumcases = Matrix{Union{Missing, Float64}}(missing,T,N)
  recovered= Matrix{Union{Missing, Float64}}(missing,T,N)
  severe= Matrix{Union{Missing,Float64}}(missing,T,N)
  name = Vector{String}(undef,N)
  day0 = minimum(skipmissing(df[!,datevar])) - Dates.Day(tstart)
  population=Vector{Float64}(undef, N)
  for (i, id) in enumerate(ids)
    dfi = df[df[!,idvar].==id,:]
    @assert length(unique(dfi[!,namevar]))==1
    name[i] = unique(dfi[!,namevar])[1]
    population[i] = unique(skipmissing(dfi[!,popvar]))[1]
    for (t, date) in enumerate(dates)
      tvals[t,i] = (date - day0).value
      inc = dfi[!,datevar].==date
      if (sum(inc)==0)
        X[:,t,i] .= 0
        continue
      end
      @assert sum(inc)==1
      inc = findfirst(inc)
      dead[t,i] = dfi[inc,deadvar]
      recovered[t,i] = dfi[inc,recoveredvar]
      cumcases[t,i] = dfi[inc,casevar] #cases - recovered[t,i] - dead[t,i]
      severe[t,i] = dfi[inc,severevar]
      k = 1
      for pv in policyvars
        d = unique(dfi[!,pv])[1]
        if d === missing
          X[k,t,i] = 0
        else
          X[k,t,i] = dates[t] >= d
        end
        k += 1
      end
      for v in otherxvars
        val=dfi[inc, v]
        if (t>60)
          mval = X[k,t-1,i]
        else
          mval = 0
        end
        X[k,t,i] = ismissing(val) ? mval : -val/100
        k += 1
      end
    end
  end
  return(RegionsData(population, tvals, dead, recovered, cumcases, severe, X, name))
end
