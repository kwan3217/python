grammar S_expr;

/*
 * Parser Rules
 */

s_epxr              : generic+ EOF ;

generic             : '(' WORD 
                       ( string
                       | generic
                       )* ')' ;

string              : WORD
                    | QUOTEDSTRING 
                    ;

/*
 * Lexer Rules
 */


WORD                : [-A-Za-z0-9_/\\.:*~?+#]+ ;
QUOTEDSTRING        : '"' ~["]* '"';

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;


