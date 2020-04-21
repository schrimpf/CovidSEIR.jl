BIB = jmd/covid.bib
PANDOC = pandoc --bibliography $(BIB) -f markdown -t markdown_mmd-citations --metadata link-citations=true

WEAVE = julia --project=@. makeweave.jl
vpath %.jmd jmd


jmdin := $(wildcard jmd/*.jmd)
jmdmd := $(patsubst jmd/%.jmd, md/%.md,$(jmdin))
jmdout := $(patsubst jmd/%.jmd, build/%.md,$(jmdin))

docin := $(wildcard src/*.md)
docout := $(patsubst src/%.md, build/%.md,$(docin))


mdout := $(jmdout)
nbout := $(patsubst jmd/%.jmd, build/%.ipynb,$(jmdin))


jmd : $(jmdout)
nb : $(nbout)
doc: $(docout)


$(jmdmd) : md/%.md : %.jmd #makeweave.jl
	$(WEAVE) $? weave

$(mdout) : build/%.md : md/%.md
	$(PANDOC) $< -o $@

$(nbout) : build/%.ipynb : %.jmd
	$(WEAVE) $^  notebook

$(docout) : build/%.md :  src/%.md
	julia --project=@. -e "using Documenter, CovidSEIR, DocumenterMarkdown; makedocs( modules=[CovidSEIR],  format=Markdown(), clean=false,  pages=[  \"Home\" => \"index.md\" ], repo=\"https://github.com/schrimpf/CovidSEIR.jl/blob/{commit}{path}#L{line}\",  sitename=\"CovidSEIR.jl\",  authors=\"Paul Schrimpf <paul.schrimpf@gmail.com>\")"

site: site/index.html

site/index.html: $(rmdout) $(jmdout) $(docout) mkdocs.yml build/figures
	mkdocs build

all: site/index.html #$(nbout)

deploy: site/index.html
	mkdocs gh-deploy


estpref = canada china italy korea us
#estpref = china italy korea us
estnew := $(patsubst %, %_tv_$(shell date +'%Y-%m.%d').jld2, $(estpref))

$(estnew):
	julia --project=@. estimate_dhmc.jl $@

estimate: $(estnew)
