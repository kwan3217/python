grammar KicadSimple;

/*
 * Parser Rules
 */

kicad               : exportnode EOF ;

exportnode          : '(' EXPORT (versionnode | content)+ ')' ;
versionnode         : '(' VERSION (WORD | STRING) ')' ;
node                : '(' WORD content* ')' ;

content             : WORD | STRING | node ;

/*
 * Lexer Rules
 */

EXPORT              : 'export';
VERSION             : 'version';
WORD                : ([-A-Za-z0-9_/\\.:*~?+] )+ ;

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;

STRING : '"' ~[<"]* '"';

