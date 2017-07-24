all: py

#This disables implicit rules:
#https://stackoverflow.com/questions/4122831/disable-make-builtin-rules-and-variables-from-inside-the-make-file
.SUFFIXES:

GRAMMAR=NetlistFile SchematicFile PcbFile S_expr
G4=$(addsuffix .g4,$(GRAMMAR))
PY=$(G4:.g4=/Parser.py)

py: $(PY)

$(GRAMMAR) : % :
	mkdir -v $@

%/Parser.py: %.g4 %
	$(eval BASE=$(shell basename $< .g4))
	touch $(BASE)/__init__.py
	antlr4 -o $(BASE) -visitor -no-listener -Dlanguage=Python3 $<

clean:
	$(RM) -rf $(GRAMMAR)
	$(RM) -rf __pycache__

#(Failed) Attempt to patch NetlistFile/NetlistFileParser.py into NetlistFile/Parser.py
#	for i in $(BASE)/$(BASE)*.py; do \
#		cat $$i \
#		| sed -e 's/\bLexer\b/antlr4.Lexer/g' \
#		| sed -e 's/\bParser\b/antlr4.Parser/g' \
#		| sed -e '/from antlr4 import/a import antlr4' \
#		| sed -e 's/$(BASE)*//g' \
#		> $(BASE)/$${i#$(BASE)/$(BASE)} ;\
#	done
#		| sed -e 's/\bATNDeserializer\b/antlr4.ATNDeserializer/' \
#		| sed -e 's/\bDFA\b/antlr4.DFA/' \
#		| sed -e 's/\bToken\b/antlr4.Token/' \
#		| sed -e 's/\bTokenStream\b/antlr4.TokenStream/' \
#		| sed -e 's/\bParserRuleContext\b/antlr4.ParserRuleContext/' \
#		| sed -e 's/\bPredictionContextCache\b/antlr4.PredictionContextCache/' \

