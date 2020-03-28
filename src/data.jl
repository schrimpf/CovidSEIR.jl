"""
    covidjhudata()

Downloads most recent JHU CSSE data on covid cases, deaths, and recoveries.

Returns a DataFrame
"""
function covidjhudata()

  url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
  ids = Symbol.(["Province/State","Country/Region", "Lat","Long"])
  dfs=map(name->begin
          res = HTTP.get(replace(url, "confirmed" => name))
          df = CSV.read(res.body)
          df = DataFrames.stack(df, DataFrames.Not(ids), variable_name=:Date, value_name = Symbol(name))
          end, ["confirmed", "deaths", "recovered"])

  df = DataFrames.join(dfs..., on=[ids..., :Date], kind=:outer)
  df[!,:Date] .= Dates.Date.(String.(df[!,:Date]), "mm/dd/yy") .+ Dates.Year(2000)
  for col in names(df)
    df[!,col]=replace(df[!,col], NaN=>missing)
  end
  DataFrames.rename!(df, [Symbol("Province/State") => :Province,
                          Symbol("Country/Region") => :Country])

  # merge country populations from world bank
  ambig = []
  df.iso2c = Vector{Union{Missing, String}}(undef, DataFrames.nrow(df))
  df.iso2c .= missing
  for c in unique(df.Country)
    res = WorldBankData.search_wdi("countries","name",Regex(c))
    if DataFrames.nrow(res) != 1
      push!(ambig, (c, res))
    else
      df[df.Country.==c, :iso2c] .= res.iso2c[1]
    end
  end
  # hand checked corrections
  df[df.Country.=="China", :iso2c] .= "CN"
  df[(df.Country.=="China") .& (df.Province.=="Hong Kong"), :iso2c] .= "HK"
  df[df.Country.=="Congo (Brazzaville)", :iso2c] .= "CG"
  df[df.Country.=="Congo (Kinshasa)", :iso2c] .= "CD"
  df[df.Country.=="Czechia", :iso2c] .= "CZ"
  df[df.Country.=="Guinea", :iso2c] .= "GN"
  df[df.Country.=="Korea, South", :iso2c] .= "KR"
  df[df.Country.=="Kyrgyzstan", :iso2c] .= "KG"
  df[df.Country.=="Niger", :iso2c] .= "CZ"
  df[df.Country.=="Nigeria", :iso2c] .= "NG"
  df[df.Country.=="Saint Lucia", :iso2c] .= "LC"
  df[df.Country.=="Saint Vincent and the Grenadines", :iso2c] .= "VC"
  df[df.Country.=="Slovakia", :iso2c] .= "SK"
  df[df.Country.=="South Africa", :iso2c] .= "ZA"
  df[df.Country.=="Sudan", :iso2c] .= "SD"
  df[df.Country.=="Dominica", :iso2c] .= "DM"

  csvfile=normpath(joinpath(dirname(Base.find_package("CovidSEIR")),"..","data","cpop.csv"))
  redownload=false
  if (redownload || !isfile(csvfile))
    tmp = DataFrames.DataFrame(iso2c=unique(skipmissing(df.iso2c)))
    tmp.cpop = Vector{Union{Missing,Float64}}(undef, DataFrames.nrow(tmp))
    tmp.cpop .= missing
    for c in tmp.iso2c
      println(c)
      tmp.cpop[c.==tmp.iso2c] .=
        try
          wb = WorldBankData.wdi("SP.POP.TOTL",c, 2010)
          sort!(wb,[:year])
          r = DataFrames.nrow(wb)
          while ismissing(wb.SP_POP_TOTL[r]) && r>1
            r -= 1
          end
          wb.SP_POP_TOTL[r]
        catch
          missing
        end
    end
    CSV.write(csvfile, tmp)
  end
  tmp = CSV.read(csvfile)
  df = join(df, tmp, on=[:iso2c], kind=:left)

  # Add province populations for Canada
  zipfile=normpath(joinpath(dirname(Base.find_package("CovidSEIR")),"..","data","17100009-eng.zip"))
  if !isfile(zipfile)
    download("https://www150.statcan.gc.ca/n1/tbl/csv/17100009-eng.zip", zipfile)
  end
  csvfile=normpath(joinpath(dirname(Base.find_package("CovidSEIR")),"..","data","17100009.csv"))
  if !isfile(csvfile)
    datadir = normpath(joinpath(dirname(Base.find_package("CovidSEIR")),"..","data"))
    run(`unzip $zipfile -d $datadir`)
  end
  cnd=CSV.read(csvfile)
  lastdate=maximum(cnd.REF_DATE)
  cnd=cnd[cnd.REF_DATE.==lastdate,[:GEO,:VALUE]]
  DataFrames.rename!(cnd, [:GEO => :Province, :VALUE => :ppop])
  df = join(df, cnd, on=[:Province], kind=:left)

  # province populations for China
  chinapop=Dict("Guangdong" => 111690000, "Shandong" => 100060000, "Henan" => 95590000,
               "Sichuan" => 83020000, "Jiangsu" => 80290000, "Hebei" => 75200000 ,
               "Hunan" => 68600000 , "Anhui" => 62550000 ,  "Hubei" => 59020000 ,
               "Zhejiang" => 56570000, "Guangxi" => 48850000, "Yunnan" =>  48010000 ,
               "Jiangxi" => 46220000 , "Liaoning" => 43690000 , "Fujian" => 39110000,
               "Shaanxi" => 38350000, "Heilongjiang" => 37890000,"Shanxi" => 36820000 ,
               "Guizhou" => 35550000,  "Chongqing" => 30750000, "Jilin" => 27170000,
               "Gansu" => 26260000 , "Inner Mongolia" => 25290000, "Xinjiang" => 24450000,
               "Shanghai" => 24180000, "Taiwan " => 23562318, "Beijing" => 21710000 ,
               "Tianjin" => 15570000 ,  "Hainan" => 9170000 , "Hong Kong Hong Kong" => 7335384 ,
               "Ningxia" => 6820000 , "Qinghai" => 5980000 ,  "Tibet" => 3370000 ,
               "Macau" => 644900  )
  for k in keys(chinapop)
    df.ppop[.!ismissing.(df.Province) .& (df.Province.==k)] .= chinapop[k]
  end
  cp = maximum(df.cpop[df.Country.=="China"])
  hk = (df.Country.=="China") .& (df.Province.=="Hong Kong")
  hkp = unique(df.cpop[hk])
  df.cpop[hk] .= cp
  df.ppop[hk] .= hkp

  # province populations for Australia
  auspop = Dict("New South Wales" => 7317500, "Victoria" => 5640900,
                "Queensland" => 4599400, "Western Australia" =>2366900,
                "South Australia" => 1659800, "Tasmania" =>511000,
                "Australian Capital Territory" => 366900,
                "Northern Territory" => 231200)
  for k in keys(auspop)
    df.ppop[.!ismissing.(df.Province) .& (df.Province.==k)] .= auspop[k]
  end

  return(df)

end

