grammar KicadSimple;

/*
 * Parser Rules
 */

kicad               : node+ EOF ;

node                : '(' WORD WHITESPACE content ')' ;

content             : TEXT | STRING | node+ ;

/*
 * Lexer Rules
 */

STRING : '"' ~[<"]* '"'

fragment LOWERCASE  : [a-z] ;
fragment UPPERCASE  : [A-Z] ;

WORD                : (LOWERCASE | UPPERCASE | '_')+ ;

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ ;

TEXT                : ('['|'(') ~[\])]+ (']'|')');
