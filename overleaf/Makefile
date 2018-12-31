RMD_FILE = ../Rmd/2018-12-16_scatterplots-13local-147nonlocal.Rmd

# Tools
LATEXMK = latexmk
RM = rm -f

# Project-specific settings
#DOCNAME = main

# Targets
all: pdf
pdf: main.pdf



# Rules
main.pdf: 
	pdflatex main; biber main; pdflatex main

mostlyclean:
	$(LATEXMK) -silent -c
	$(RM) *.bbl

clean: mostlyclean
	$(LATEXMK) -silent -C
	$(RM) *.run.xml *.synctex.gz
	$(RM) *.d

.PHONY: all clean mostlyclean pdf


