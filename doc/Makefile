css_style = jbldoc.css
# Make sure to have computer modern fonts and lyx fonts installed, otherwise
#   you are in for a world of pain...
# Besides the obvious dependencies, which make the compilation fail, check that
#   you have something like lyx-fonts and fonts-lmodern installed.
# See https://github.com/asciidoctor/asciidoctor-mathematical for details.
opts_pdf = -a toc="left" -a toclevels="3" -r asciidoctor-mathematical -a mathematical-format=pdf
opts_html = -a stylesheet=$(css_style) -a toc="left" -a toclevels="3"
opts_epub = -a toc="left" -a toclevels="3"

EXT = html

prog := $(wildcard prog_*.adoc)
prog_ext := $(prog:.adoc=.$(EXT))

lib =  estampes_units.adoc estampes_lib.adoc
lib_ext := $(lib:.adoc=.$(EXT))

api = estampes_api.adoc
api_ext := $(api:.adoc=.$(EXT))

main = estampes.adoc
main_ext := $(main:.adoc=.$(EXT))

full = $(main) $(prog) $(lib) $(api)
full_ext := $(full:.adoc=.$(EXT))
full_pdf := $(full:.adoc=.pdf)
full_html := $(full:.adoc=.html)
full_epub := $(full:.adoc=.epub)


all: $(full_ext)

full: $(full_pdf) $(full_html)

estampes: $(main_ext)

progs: $(prog_ext)

api: $(api_ext)

lib: $(lib_ext)

autodoc :
	python3 $$(which sphinx-apidoc) -f -o source/ ../../estampes/ ../../estampes/setup.py
	python3 $$(which sphinx-build) -b html source/ build/

%.pdf : %.adoc
	asciidoctor-pdf $< $(opts_pdf) -o build/pdf/$@
	rm stem-*

%.html : %.adoc
	asciidoctor $< $(opts_html) -o build/html/$@

%.epub : %.adoc
	asciidoctor-epub3 $< $(opts_epub) -o build/epub/$@

clean:
	rm -fr build/

autoclean :
	rm source/estampes*.rst source/modules.rst source/test.rst
	rm -r build/*
