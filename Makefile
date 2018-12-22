RMD_FILE = Rmd/2018-12-16_scatterplots-13local-147nonlocal.Rmd

# Tools
LATEXMK = latexmk
RM = rm -f

# Project-specific settings
DOCNAME = overleaf/main

# Targets
all: doc
doc: pdf
pdf: overleaf/main.pdf
fig: Rmd/%.eps
data: data/%.rds data-to-ignore/%.Rdata


# Rules
%.pdf: %.tex fig
	cd overleaf; $(LATEXMK) -pdfps -bibtex main.tex
mostlyclean:
	cd overleaf; $(LATEXMK) -silent -c
	$(RM) overleaf/*.bbl

clean: mostlyclean
	cd overleaf; $(LATEXMK) -silent -C
	$(RM) overleaf/*.run.xml overleaf/*.synctex.gz
	$(RM) overleaf/*.d

.PHONY: all clean doc mostlyclean pdf

# Include auto-generated dependencies
#-include *.d



	
Rmd/%.eps: $(RMD_FILE)
	Rscript -e "rmarkdown::render('$<')"
  
