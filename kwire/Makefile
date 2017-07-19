all: py

#This disables implicit rules:
#https://stackoverflow.com/questions/4122831/disable-make-builtin-rules-and-variables-from-inside-the-make-file
.SUFFIXES:

GRAMMAR=Netlist Schematic
G4=$(addsuffix .g4,$(GRAMMAR))
PY=$(G4:.g4=/Parser.py)

py: $(PY)

$(GRAMMAR) : % :
	mkdir -v $@

%/Parser.py: %.g4 %
	touch `basename $< .g4`/__init__.py
	antlr4 -o `basename $< .g4` -visitor -no-listener -Dlanguage=Python3 $<

clean:
	$(RM) -rf $(GRAMMAR)
	$(RM) -rf __pycache__

