grammar PcbFile;

/*
 * Parser Rules
 */

pcb                 : kicad_pcb EOF ;

kicad_pcb           : KICAD_PCB
                       ( netdef
                       | module
                       | gr_arc
                       | gr_text
                       | gr_line
                       | segment
                       | via
                       | zone
                       | generic
                       )+ ')' ;
module              : MODULE string 
                       (layer
                       |at
                       |fp_text
                       |fp_line
                       |pad
                       |generic
                       )+ ')';
gr_arc              : GR_ARC
                       (start
                       |end
                       |angle
                       |layer
                       |width
                       |generic
                       )+ ')';
gr_text             : GR_TEXT string
                       (at
                       |layer
                       |effects
                       |generic
                       )+ ')';
gr_line             : GR_LINE
                       (start
                       |end
                       |layer
                       |width
                       |generic
                       )+ ')';
segment             : SEGMENT
                       (start
                       |end
                       |layer
                       |width
                       |netref
                       |generic
                       )+ ')';
via                 : VIA
                       (at
                       |size1
                       |drill
                       |layers
                       |netref
                       |generic
                       )+ ')';
zone                : ZONE
                       (netref
                       |layer
                       |filled_polygon
                       |generic
                       )+ ')';
fp_text             : FP_TEXT string string
                       (at
                       |layer
                       |effects
                       |generic
                       )+ ')';
fp_line             : FP_LINE
                       (start
                       |end
                       |layer
                       |width
                       |generic
                       )+ ')';
pad                 : PAD string string string
                       (at
                       |layers
                       |size2
                       |netdef
                       |drill
                       )+ ')';
effects             : EFFECTS 
                       (font
                       |generic
                       )+ ')';
layers              : LAYERS string string+ ')' ;
filled_polygon      : FILLED_POLYGON pts ')' ;
pts                 : PTS xy+ ')' ;                     
font                : FONT
                       (size2
                       |thickness
                       )+ ')';
/* Generic node - capable of holding subnodes. Used for don't-care nodes */
generic             : '(' WORD 
                       ( string
                       | generic
                       )* ')' ;

/* Collection nodes. These are collections of the same type of node. These 
   collections will have pluralized names and will be collections of nodes
   with the corresponding singular name. */
/*
components      : COMPONENTS comp+    ')' ;
fields          : FIELDS     field+   ')';
footprints      : FOOTPRINTS fp+      ')';
libraries       : LIBRARIES  library+ ')' ;
libparts        : LIBPARTS   libpart+ ')' ;
nets            : NETS       net+     ')' ;
pins            : PINS       pin+     ')' ;
aliases         : ALIASES    alias+   ')' ;
*/

/* Boring nodes - each only is allowed to have a single string value and no subnodes */
size1           : SIZE      string ')' ;
drill           : DRILL     string ')' ;
angle           : ANGLE     string ')' ;
netref          : NET       string ')' ;
width           : WIDTH     string ')' ;
layer           : LAYER     string ')' ;
thickness       : THICKNESS string ')' ;

/* Double boring nodes - twice as boring, must have two strings */
netdef          : NET    string string ')' ;
start           : START  string string ')' ;
end             : END    string string ')' ;
xy              : XY     string string ')' ;
size2           : SIZE   string string ')' ;

/* Triply boring nodes - must have three strings */
at              : AT    string string string')' ;

string              : WORD
                    | QUOTEDSTRING 
                    ;

/*
 * Lexer Rules
 */

/* Node names. Includes the parenthesis at the beginning so that we don't get
   confused when a part is named alias and we get (name alias) */
KICAD_PCB           : '(kicad_pcb';
MODULE              : '(module';
GR_ARC              : '(gr_arc';
GR_LINE             : '(gr_line';
GR_TEXT             : '(gr_text';
SEGMENT             : '(segment';
VIA                 : '(via';
ZONE                : '(zone';
FP_TEXT             : '(fp_text';
FP_LINE             : '(fp_line';
PAD                 : '(pad';
EFFECTS             : '(effects';
LAYERS              : '(layers';
FILLED_POLYGON      : '(filled_polygon';
PTS                 : '(pts';
FONT                : '(font';
SIZE                : '(size';
DRILL               : '(drill';
ANGLE               : '(angle';
NET                 : '(net';
WIDTH               : '(width';
LAYER               : '(layer';
THICKNESS           : '(thickness';
START               : '(start';
END                 : '(end';
XY                  : '(xy';
AT                  : '(at';

WORD                : [-A-Za-z0-9_/\\.:*~?+#]+ ;
QUOTEDSTRING        : '"' ~["]* '"';

WHITESPACE          : (' ' | '\t' | '\r' | '\n')+ -> skip;


