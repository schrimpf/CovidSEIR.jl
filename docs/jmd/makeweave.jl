using Pkg
Pkg.activate("..")
using Weave

runweave= true
runnotebook= false

src="covid.jmd"

if runweave
  println("weaving markdown for $src")
  using Weave
  weave(src,out_path=".",
        cache=:user, cache_path="weavecache",
        doctype="github", mod=Main,
        args=Dict("md" => true))
  pdin = "covid.md"
  pdout = "../src/covid.md"
  BIB = "covid.bib"
  run(`pandoc --bibliography $(BIB) -f markdown -t markdown_mmd-citations --metadata link-citations=true $pdin -o $pdout`)
end

if runnotebook
  println("weaving notebook for $src")
  using Weave
  notebook(src, out_path=joinpath(pwd(),"../src"), nbconvert_options="--allow-errors")
end

