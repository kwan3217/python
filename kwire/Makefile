test_%: %Parser.class
	$(eval BASENAME=$(shell echo $@ | sed 's/test_//'))
	$(eval FIRSTRULE=$(shell grep : $(BASENAME).g4 | head -n 1 | sed 's/\s*:.*$///'))
	grun $(BASENAME) $(FIRSTRULE) -gui < $(TESTDATA)

all: class py

GRAMMAR=Netlist.g4 Schematic.g4
JAVA=$(GRAMMAR:.g4=Parser.java)
CLASS=$(GRAMMAR:.g4=Parser.class)
PY=$(GRAMMAR:.g4=Parser.py)
TESTDATA=../../kicad/Tomo/Tomo.net

class: $(CLASS) java
java: $(JAVA)
py: $(PY)

%Parser.java: %.g4
	antlr4 $<

%Parser.class: %Parser.java
	javac $(^:Parser.java=*.java)

%Parser.py: %.g4
	antlr4 -visitor -no-listener -Dlanguage=Python3 $<

clean: 
	$(RM) $(JAVA:Parser.java=Parser*.java)
	$(RM) $(JAVA:Parser.java=Lexer*.java)
	$(RM) $(JAVA:Parser.java=BaseListener*.java)
	$(RM) $(JAVA:Parser.java=Listener*.java)
	$(RM) $(CLASS:Parser.class=Parser*.class)
	$(RM) $(CLASS:Parser.class=Lexer*.class)
	$(RM) $(CLASS:Parser.class=BaseListener*.class)
	$(RM) $(CLASS:Parser.class=Listener*.class)
	$(RM) $(GRAMMAR:.g4=.tokens)
	$(RM) $(GRAMMAR:.g4=Lexer.tokens)
	$(RM) $(GRAMMAR:.g4=Parser.py)
	$(RM) $(GRAMMAR:.g4=Lexer.py)
	$(RM) $(GRAMMAR:.g4=Listener.py)
	$(RM) $(GRAMMAR:.g4=Visitor.py)
	$(RM) -rf __pycache__
