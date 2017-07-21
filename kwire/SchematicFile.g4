grammar Schematic;

/*
 * Parser Rules
 */

schematic           : header comps EOF ;

string              : QUOTEDSTRING | WORD ;

header              : signature libss eelayer descr;

signature           : SIGNATURE INT ;

libss               : (libs)* ;

libs                : 'LIBS:' string ;

eelayer             : 'EELAYER' INT INT 'EELAYER' 'END';

descr               : '$Descr' 
                    (sheetsize 
                    |encoding
                    |sheet
                    |title
                    |date
                    |rev
                    |descrcomp
                    |comments
                    )+ '$EndDescr';
sheetsize           : string INT INT;
encoding            : 'encoding' string;
sheet               : 'Sheet' INT INT;
title               : 'Title' string;
date                : 'Date'  string;
rev                 : 'Rev'   string;
descrcomp           : 'Comp'  string;
comments            : 'Comment1' string 'Comment2' string 'Comment3' string 'Comment4' string ;

comps               : (comp)*;
comp                : '$Comp' 
                    (compID
                    |compU
                    |compPos
                    |compAlternate
                    |compField
                    )+ '$EndComp' ;
compID              : 'L' string string ;
compU               : 'U' INT INT TIMESTAMP ;
compPos             : 'P' INT INT ;
compAlternate       : 'AR' 'Path=' string 'Ref=' string 'Part=' string ;
compField           : 'F' INT string string INT INT INT INT string string compFieldMatrix? ;
compFieldMatrix     : '1' INT INT INT INT INT INT ;

/*
 * Lexer Rules
 */

/* Node names. Includes the parenthesis at the beginning so that we don't get
   confused when a part is named alias and we get (name alias) */

SIGNATURE           : 'EESchema Schematic File Version' ;
TIMESTAMP           : [0-9A-F][0-9A-F][0-9A-F][0-9A-F][0-9A-F][0-9A-F][0-9A-F][0-9A-F] ;
INT                 : ('+'|'-')?[0-9]+ ;

QUOTEDSTRING        : ('"' ~["]* '"') ;
WORD                : [-A-Za-z0-9_/\\.*~?+]+ ;

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;


