#
# Makefile
#
# $Id: Makefile 1119 2008-11-04 21:23:52Z mark $
#
# Global project Makefile 
#
#
# This makefile was based on
# http://209.85.135.104/search?q=cache:-fU4UEn36AgJ:kriemhild.uft.uni-bremen.de/viewcvs/Makefile.rnoweb+makefile+sweave+latex&hl=en&ct=clnk&cd=10&client=safari
# (The original site found on Google is no longer available and I couldn't track a new version)
# This is credited to Johannes Ranke, based on work by Nicholas Lewin-Koh and Rouben Rostmaian
# Further advice was provided by Deepayan Sarkar via the R-help mailing list.
#


# The master LaTeX document is thesis.tex, and \include{}'s the required chapters
# Other master documents are acnr_master.tex (which includes acnr_paper.tex), genetic_master and drpla_master
# This approach ensures integrity of the Sweave files for papers
MASTER = thesis.pdf acnr_master.pdf genetic_master.pdf genetic_suppl.pdf drpla_master.pdf icars_master.pdf clinical_master.pdf drplarv_master.pdf

# the master document depends on all of the tex files
RNWFILES = $(wildcard *.Rnw)
TEXFILES = $(wildcard *.tex)
DEPENDS = $(patsubst %.Rnw,%.tex,$(RNWFILES)) $(TEXFILES) $(RFILES) 


RERUN = "(There were undefined references|Rerun to get (citations|cross-references|the bars) (correct|right)|Table widths have changed. Rerun LaTeX.|Linenumber reference failed)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined"


# density to use for conversion of vector graphics into bitmaps (rasterisation density)
DENSITY = 150
VECTOR_OMIT = "thesis.pdf|*_master.pdf"	# don't convert certain specific PDFs


all: $(MASTER)

view: $(MASTER)
	@open $(MASTER)
	@osascript -e 'tell application "System Events"' -e 'tell process "Preview"' -e 'tell menu bar 1' -e 'tell menu bar item "File"' -e 'tell menu 1' -e 'click menu item "Close"' -e 'end tell' -e 'end tell' -e 'end tell' -e 'end tell' -e 'end tell'
	@open $(MASTER)

$(MASTER): $(DEPENDS)

%.tex: %.Rnw ataxia.R
	R --arch i386 CMD SWEAVE '$<'

%.pdf: %.tex
	@pdflatex $<
	@makeindex -s $*.ist -t $*.alg -o $*.acr $*.acn
	@egrep -c $(RERUNBIB) $*.log && (bibtex $*;pdflatex $<); true
	@egrep $(RERUN) $*.log && (pdflatex $<) ; true
	@egrep $(RERUN) $*.log && (pdflatex $<) ; true

%.txt: %.pdf
	@pdftotext -nopgbrk $<

clean:
	@rm -f *.aux *.log *.bbl *.blg *.brf *.cb *.ind *.idx *.ilg  \
          *.inx *.ps *.dvi *.toc *.out *.lot *~ *.lof *.alg *.acn *.acr *.glo *.ist *.ttt *.fff \
          *.eps *.pdf *.svn *.dot *.lol
	@rm -f $(patsubst %.Rnw,%.tex,$(RNWFILES))
	@rm -rf dynamic_jpg dynamic_png

# some helper rules to generate JPG or PNG versions of all of the graphics
jpeg:	all
	@mkdir -p dynamic_jpg
	@for i in `ls *.pdf | egrep -v $(VECTOR_OMIT)`; do echo "$$i"; convert -density $(DENSITY) "$$i" "dynamic_jpg/$$i.jpg"; done

png:	all
	@mkdir -p dynamic_png
	@for i in `ls *.pdf | egrep -v $(VECTOR_OMIT)`; do echo "$$i"; convert -density $(DENSITY) "$$i" "dynamic_png/$$i.png"; done
