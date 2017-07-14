test: all
	grun $(GRAMMAR:.g4=) `grep : $(GRAMMAR) | head -n 1 | sed 's/\s*:.*$///'` -gui < $(TESTDATA)

all: class py

GRAMMAR=Schematic.g4
JAVA=$(GRAMMAR:.g4=.java)
CLASS=$(GRAMMAR:.g4=.class)
PY=$(GRAMMAR:.g4=Listener.py)
TESTDATA=1bitadder.sch

class: $(CLASS) java
java: $(JAVA)
py: $(PY)

%.java: %.g4
	antlr4 $<

%.class: %.java
	javac $(^:.java=*.java)

%Listener.py: %.g4
	antlr4 -Dlanguage=Python3 $<

clean: 
	$(RM) $(JAVA:.java=*.java)
	$(RM) $(CLASS:.class=*.class)
	$(RM) $(GRAMMAR:.g4=.tokens)
	$(RM) $(GRAMMAR:.g4=Lexer.tokens)
	$(RM) $(GRAMMAR:.g4=Parser.py)
	$(RM) $(GRAMMAR:.g4=Lexer.py)
	$(RM) $(GRAMMAR:.g4=Listener.py)

