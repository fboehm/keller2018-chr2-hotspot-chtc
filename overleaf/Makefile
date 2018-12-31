RMD_FILE = Rmd/2018-12-16_scatterplots-13local-147nonlocal.Rmd

# Tools
LATEXMK = latexmk
RM = rm -f

# Project-specific settings
DOCNAME = overleaf/main

# Targets
all: pdf
pdf: overleaf/main.pdf



# Rules
overleaf/main.pdf: 
	cd overleaf; $(LATEXMK) -pdfps -bibtex main.tex

mostlyclean:
	cd overleaf; $(LATEXMK) -silent -c
	$(RM) overleaf/*.bbl

clean: mostlyclean
	cd overleaf; $(LATEXMK) -silent -C
	$(RM) overleaf/*.run.xml overleaf/*.synctex.gz
	$(RM) overleaf/*.d

.PHONY: all clean mostlyclean pdf


