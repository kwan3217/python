grammar Kicad;

/*
 * Parser Rules
 */

kicad               : exportnode EOF ;

exportnode          : '(' EXPORT genericnode* componentsnode genericnode* ')' ;

componentsnode      : '(' COMPONENTS compnode* ')' ;

compnode            : '(' COMP compinternals* ')' ;

refnode             : '(' REF (WORD | STRING) ')' ;

valuenode           : '(' VALUE (WORD | STRING) ')' ;

libsourcenode       : '(' LIBSOURCE '(' LIB (WORD|STRING) ')' '(' PART (WORD|STRING) ')' genericnode* ')' ;

compinternals       : refnode 
                    | valuenode
                    | libsourcenode
                    | genericnode
                    ;

node                : genericnode 
                    | componentsnode 
                    | compinternals
                    ;

genericnode         : '(' WORD ')'
                    | '(' WORD content ')' ;

content             : (WORD | STRING | node )+ ;


/*
 * Lexer Rules
 */

EXPORT              : 'export';
COMPONENTS          : 'components';
COMP                : 'comp';
REF                 : 'ref';
VALUE               : 'value';
LIBSOURCE           : 'libsource';
LIB                 : 'lib';
PART                : 'part';

WORD                : ([-A-Za-z0-9_/.:*~?+] )+ ;

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;

STRING              : '"' ~["]* '"';




