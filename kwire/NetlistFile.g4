grammar NetlistFile;

/*
 * Parser Rules
 */

netlist             : export EOF ;

export              : EXPORT 
                       (version 
                       | design
                       | components
                       | libparts
                       | libraries
                       | nets
                       )+ ')' ;
design              : DESIGN 
                       (source
                       |date
                       |tool
                       |sheet
                       )+ ')';
sheet               : SHEET 
                       (number
                       |name
                       |tstamps
                       |title_block
                       )+ ')';
title_block         : TITLE_BLOCK
                       (title
                       |company
                       |rev
                       |date
                       |source
                       |comment
                       )+ ')';
comment             : COMMENT
                       (number
                       |value
                       )+ ')';
comp                : COMP
                       (ref
                       |value
                       |footprint
                       |fields
                       |libsource
                       |sheetpath
                       |tstamp
                       )+ ')';
libpart             : LIBPART
                       (lib
                       |part
                       |aliases
                       |description
                       |docs
                       |fields
                       |pins
                       |footprints
                       )+ ')';
pin                 : PIN
                       (num
                       |name
                       |type_
                       )+ ')';
net                 : NET
                       (code
                       |name
                       |node
                       )+ ')';
node                : NODE
                       (ref
                       |pin__innode
                       )+ ')';
field               : FIELD name string ')';
library             : LIBRARY 
                       (logical 
                       |uri
                       )+ ')';        
libsource           : LIBSOURCE lib part ')';        
sheetpath           : SHEETPATH names tstamps ')';        

/* Generic node - capable of holding subnodes */
/*
generic             : '(' WORD 
                       ( string
                       | generic
                       )* ')' ;
*/

/* Collection nodes. These are collections of the same type of node. These 
   collections will have pluralized names and will be collections of nodes
   with the corresponding singular name. */
components      : COMPONENTS comp+    ')' ;
fields          : FIELDS     field+   ')';
footprints      : FOOTPRINTS fp+      ')';
libraries       : LIBRARIES  library+ ')' ;
libparts        : LIBPARTS   libpart+ ')' ;
nets            : NETS       net+     ')' ;
pins            : PINS       pin+     ')' ;
aliases         : ALIASES    alias+   ')' ;

/* Boring nodes - each only is allowed to have a single string value and no subnodes */
alias           : ALIAS       string ')' ;
code            : CODE        string ')' ;
description     : DESCRIPTION string ')' ;
docs            : DOCS        string ')' ;
footprint       : FOOTPRINT   string ')' ;
fp              : FP          string ')' ;
lib             : LIB         string ')' ;
logical         : LOGICAL     string ')' ;
name            : NAME        string ')' ;
names           : NAMES       string ')' ;
num             : NUM         string ')' ;
number          : NUMBER      string ')' ;
part            : PART        string ')' ;
pin__innode     : PIN         string ')' ;
ref             : REF         string ')' ;
source          : SOURCE      string ')' ;
tool            : TOOL        string ')' ;
tstamp          : TSTAMP      string ')' ;
tstamps         : TSTAMPS     string ')' ;
type_           : TYPE        string ')' ;
uri             : URI         string ')' ;
value           : VALUE       string ')' ;
version         : VERSION     string ')' ;

/* Optional boring nodes - might not even have a string */
company         : COMPANY     string? ')' ;
date            : DATE        string? ')' ;
rev             : REV         string? ')' ;
title           : TITLE       string? ')' ;

string              : WORD
                    | QUOTEDSTRING 
                    ;

/*
 * Lexer Rules
 */

/* Node names. Includes the parenthesis at the beginning so that we don't get
   confused when a part is named alias and we get (name alias) */
ALIAS               : '(alias';
ALIASES             : '(aliases';
CODE                : '(code';
COMMENT             : '(comment';
COMP                : '(comp';
COMPANY             : '(company';
COMPONENTS          : '(components';
DATE                : '(date';
DESCRIPTION         : '(description';
DESIGN              : '(design';
DOCS                : '(docs';
EXPORT              : '(export';
FIELD               : '(field';
FIELDS              : '(fields';
FOOTPRINT           : '(footprint';
FP                  : '(fp';
FOOTPRINTS          : '(footprints';
LIB                 : '(lib';
LIBPART             : '(libpart';
LIBPARTS            : '(libparts';
LIBRARIES           : '(libraries';
LIBRARY             : '(library';
LIBSOURCE           : '(libsource';
LOGICAL             : '(logical';
NAME                : '(name';
NAMES               : '(names';
NET                 : '(net';
NETS                : '(nets';
NODE                : '(node';
NUM                 : '(num';
NUMBER              : '(number';
PART                : '(part';
PIN                 : '(pin';
PINS                : '(pins';
REF                 : '(ref';
REV                 : '(rev';
SHEET               : '(sheet';
SHEETPATH           : '(sheetpath';
SOURCE              : '(source';
TITLE_BLOCK         : '(title_block';
TITLE               : '(title';
TOOL                : '(tool';
TSTAMP              : '(tstamp';
TSTAMPS             : '(tstamps';
TYPE                : '(type';
URI                 : '(uri';
VALUE               : '(value';
VERSION             : '(version';

WORD                : [-A-Za-z0-9_/\\.:*~?+#]+ ;
QUOTEDSTRING        : '"' ~["]* '"';

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;


