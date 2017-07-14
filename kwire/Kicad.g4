grammar Kicad;

/*
 * Parser Rules
 */

kicad               : exportnode EOF ;

exportnode          : '(' EXPORT 
                       (versionnode 
                       | designnode
                       | node)+ ')' ;
versionnode         : '(' VERSION string ')' ;
designnode          : '(' DESIGN 
                       (sourcenode
                       |node
                       )+ ')';
sourcenode          : '(' SOURCE string ')' ;
datenode            : '(' DATE string ')' ;
node                : '(' WORD 
                       ( string
                       | node
                       )* ')' ;
string              : WORD
                    | QUOTEDSTRING 
                    ;

/*
 * Lexer Rules
 */

EXPORT              : 'export';
VERSION             : 'version';
SOURCE              : 'source';
DESIGN              : 'design';
DATE                : 'date';
WORD                : [-A-Za-z0-9_/\\.:*~?+]+ ;
QUOTEDSTRING        : '"' ~["]* '"';

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;


