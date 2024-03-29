import re
import xml.etree.ElementTree as ET
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import tex
import subprocess
import os

spaceparse=re.compile(r"\s+")
position=re.compile(r"!--Position: ([A-Za-z0-9_]+) \(([0-9]+),([0-9]+)\)")

def tex_head(oufn):
    ouf=open(oufn+'.tex','w')
    ouf.write("\\documentclass{minimal}\n")
    ouf.write("\\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}\n")
    ouf.write("\\begin{document}\n")
    ouf.write("$$\n")
    return ouf

def tex_eqn(ouf,tag,text,phantom=False):
    ouf.write("\\zsavepos{"+tag+"}")
    if phantom:
        ouf.write("\phantom{")
    ouf.write("{"+text+"}")
    if phantom:
        ouf.write("}")
    ouf.write("\n")

def tex_foot(ouf,tags):
    ouf.write("$$\n")
#    for tag in tags:
#        ouf.write("\\typeout{!--Position: "+tag+" (\\zposx{"+tag+"},\\zposy{"+tag+"})}\n")
    ouf.write("\\end{document}\n")
    ouf.close()

def tex_render(oufn,clean=False):
    # Get the coordinates of the piece by itself
    subprocess.call("latex "+oufn+".tex > /dev/null",shell=True)
    # Delete the .tex since we no longer need it
    if clean: os.remove(oufn+'.tex')
    # clean up other files we don't care about
    if clean: os.remove(oufn+'.aux')
    if clean: os.remove(oufn+'.log')
    # Render the piece
    subprocess.call("dvisvgm -e -n -bmin --keep "+oufn+".dvi > /dev/null 2>&1",shell=True)
    # Delete the DVI now that it is rendered
    if clean: os.remove(oufn+'.dvi')
    # Slurp the SVG for this equation
    with open(oufn+'.svg',"r") as inf:
        content=inf.readlines()
    # Delete the SVG
    if clean: os.remove(oufn+'.svg')
    return content


class svgdraw:
    def pathEnd(self,path):
        pathStart=(0,0)
        result=(0,0)
        for i,code in enumerate(path["code"]):
            if code==Path.MOVETO:
                pathStart=path["coords"][i]
                result=pathStart
            elif code==Path.CLOSEPOLY:
                result=pathStart
            else:
                result=path["coords"][i]
        return result
    def smoothPt(self,path):
        if(len(path["code"])<2):
            return self.pathEnd(path)
        if(path["code"][-2]!=Path.CURVE4):
            return self.pathEnd(path)
        return path["coords"][-2]
    def finishCommand(self,path,command,coords):
        if command=="M":
            path["code"].append(Path.MOVETO)
            path["coords"].append(coords)
        elif command=="L":
            pe=self.pathEnd(path)
            path["code"].append(Path.LINETO)
            path["coords"].append(coords)
        elif command=="H":
            pe=self.pathEnd(path)
            path["code"].append(Path.LINETO)
            path["coords"].append((coords[0],pe[1]))
        elif command=="V":
            pe=self.pathEnd(path)
            path["code"].append(Path.LINETO)
            path["coords"].append((pe[0],coords[0]))
        elif command=="C":
            path["code"].append(Path.CURVE4)
            path["code"].append(Path.CURVE4)
            path["code"].append(Path.CURVE4)
            path["coords"].append(coords[0:2])
            path["coords"].append(coords[2:4])
            path["coords"].append(coords[4:6])
        elif command=="S":
            sp=self.smoothPt(path)
            pe=self.pathEnd(path)
            path["code"].append(Path.CURVE4)
            path["code"].append(Path.CURVE4)
            path["code"].append(Path.CURVE4)
            dsp=(pe[0]-sp[0],pe[1]-sp[1])
            p1=(pe[0]+dsp[0],pe[1]+dsp[1])
            path["coords"].append(p1)
            path["coords"].append(coords[0:2])
            path["coords"].append(coords[2:4])
        elif command=="Z":
            path["code"].append(Path.CLOSEPOLY)
            path["coords"].append((0,0))
        else:
            raise Exception("Unhandled command "+command)
    def draw(self,ax,coords,fade):
        for use in self.uses:
            commands=self.paths[use["id"]]["code"]
            rcoords =self.paths[use["id"]]["coords"]
            ref=use["coords"]
            acoords=[]
            for rcoord in rcoords:
                acoords.append((rcoord[0]+ref[0]+coords[0],rcoord[1]+ref[1]+coords[1]))
            path=Path(acoords,commands)
            patch=patches.PathPatch(path,facecolor=(0,0,0,fade),lw=0)
            ax.add_patch(patch)
    def set_svg(self,svg):
        self.svg=svg
        #This parser is not for general purpose SVG files, but for those 
        #generated by dvisvgm. This uses paths defined in the defs section
        #and instantiated in a use section. The path might be disconnected
        self.paths={}
        self.uses=[]
        #Parse the XML into an ElementTree
        root=ET.fromstringlist(self.svg)
        bbox=spaceparse.split(root.attrib["viewBox"])
        self.left=float(bbox[0])
        self.top=float(bbox[1])
        self.width=float(bbox[2])
        self.height=float(bbox[3])
        
        defs=root.find('{http://www.w3.org/2000/svg}defs')
        if defs is not None:
            #Get all the paths out of the xml, which are in svg/defs
            for svgpath in defs.findall("{http://www.w3.org/2000/svg}path"):
                id=svgpath.attrib["id"]
                d=re.findall('[A-Za-z ]|(?:[-+]?[0-9]+(?:\.[0-9]+)?)',svgpath.attrib["d"])
                d=[x for x in d if x != " "]
                path={"code":[],"coords":[]}
                coords=None
                command=None
                for part in d:
                    if re.match("[A-Za-z]",part):
                        if coords is not None:
                            self.finishCommand(path,command,coords)
                        command=part
                        coords=None
                    else:
                        if coords is None:
                            coords=(float(part),)
                        else:
                            coords=coords+(float(part),)
                if command is not None:
                    self.finishCommand(path,command,coords)
                self.paths[id]=path
        g=root.find('{http://www.w3.org/2000/svg}g')
        if g is not None:
            #Get all the paths out of the xml, which are in svg/defs
            for svguse in g.findall("{http://www.w3.org/2000/svg}use"):
                use={"id":svguse.attrib["{http://www.w3.org/1999/xlink}href"][1:],
                     "coords":(float(svguse.attrib["x"]),float(svguse.attrib["y"]))}
                self.uses.append(use)
    def shift(self,delta_coords):
        dx=delta_coords[0]
        dy=delta_coords[1]
        self.left+=dx
        self.top+=dy
        for use in self.uses:
            use["coords"]=(use["coords"][0]+dx,use["coords"][1]+dy)
            
class subexpr(svgdraw):
    def __init__(self,eqnname,names,texes,name):
        self.name=name
        ouf=tex_head(eqnname+name)
        for iname,itex in zip(names,texes):
            tex_eqn(ouf,iname,itex,phantom=(name!=iname))
            if(name==iname):
              self.tex =itex
        tex_foot(ouf,[name])
        svg=tex_render(eqnname+name)
        self.set_svg(svg)

class expr(svgdraw):
    def __init__(self,eqnname,names,texes,anchor):
        """
        Construct an expression out of a number of named subexpressions
        :param eqnname:  Name to use for the whole expression, used as the root
                         of the temporary files created in the process 
        :type  eqnname:  string
        :param names:  Names to use for each subexpression
        :type  names:  list of strings
        :param texes:  TeX code for each subexpression
        :type  texes:  list of strings
        :param anchor:  Name of subexpression used as anchor. In most dances,
                        the anchor subexpression won't move, and specifying the
                        location of the expression specifies the location
                        of the anchor subexpression
        :type  anchor:  list of strings

        """
        self.subexprs={}
        for name,tex in zip(names,texes):
            self.subexprs[name]=subexpr(eqnname,names,texes,name)
        #Find the original position of the subexpression
        anchor_coords=(-self.subexprs[anchor].left,-self.subexprs[anchor].top)
        #Move the subexpressions so that the anchor subexpression is at 0,0
        for k,v in self.subexprs.items():
            self.subexprs[k].shift(anchor_coords)
        pass

def cubic(t):
    if(t<0):
        return 0
    if(t>1):
        return 1
    return 6*(t*t/2-t*t*t/3)

def linterp(x0,y0,x1,y1,x):
    t=(x-x0)/(x1-x0)
    return y0*(1-t)+y1*t

class dance_eqn:
    def draw_piece(self,ax,k,t,fade=1):
        (xa,ya,fadea,xb,yb,fadeb)=self.get_piece_choreography(k,t)
        if fadea>0:
            self.exprA.subexprs[k].draw(ax,(xa,ya),fade*fadea)
        if fadeb>0:
            self.exprB.subexprs[k].draw(ax,(xb,yb),fade*fadeb)
    def draw(self,ax,t,anchor_coords,fade=1):
        for k in self.exprA.subexprs.keys():
            self.draw_piece(ax,k,t,fade)

class dance_addsub_left_right(dance_eqn):
    def __init__(self,
                 left_before ,left_sign ,term,left_after ,
                 equals,
                 right_before,right_sign,     right_after):
        """
        Render the movement of a term from the left side to the right side.
    
        :param left_before:  Part of equation on the left side to the left of 
                             (before) the term to be moved, in TeX format
        :type  left_before:  string
        :param left_sign:    Sign of the term to be moved when it is on the left side
        :type  left_sign:    string
        :param term:         Term to be moved
        :type  term:         string
        :param left_after:   Part of equation on the left side to the right of 
                             (after) the term to be moved
        :type  left_after:   string
        :param equals:       Anchor of expression, usually equals. This will remain 
                             in the same place
        :type  equals:       string
        :param right_before: Part of equation on the right side to the left of 
                             (before) the term to be moved, in TeX format
        :type  right_before: string
        :param right_sign:   Sign of the term to be moved when it is on the right side
        :type  right_sign:   string
        :param right_after:  Part of equation on the right side to the right of 
                             (after) the term to be moved
        :type  right_after:  string
        """
        # Get the spacing of the equation before the move
        self.exprA=expr("exprA",
                        ["left_before",     "sign","term","left_after","equals","right_before","right_after"],
                        [ left_before , left_sign , term , left_after , equals , right_before , right_after ],"equals")
        # Get the spacing of the equation after the move
        self.exprB=expr("exprB",
                        ["left_before","left_after","equals","right_before",      "sign","term","right_after"],
                        [ left_before , left_after , equals , right_before , right_sign , term , right_after ],"equals")
        
    def get_piece_choreography(self,k,t):
        firstStageEnd=0.2
        secondStageEnd=0.8
        xa=0
        ya=0
        fadea=1
        xb=0
        yb=0
        fadeb=1
        if(t<firstStageEnd):
            this_t=linterp(0,0,firstStageEnd,1,t)
            # No second part at all yet
            xb=0
            yb=0
            fadeb=0
            # first part completely visible
            fadea=1
            if(k=='term') or (k=='sign'):
                #move the first term and sign straight up
                xa=0
                ya=cubic(this_t)*-25
            else:
                xa=0
                ya=0
        elif(t<secondStageEnd):
            #No second part by default for most parts
            xb=0
            yb=0
            fadeb=0
            # first part completely visible
            fadea=1
            # Move parts by default
            this_t=linterp(firstStageEnd,0,secondStageEnd,1,t)
            mx=self.exprB.subexprs[k].left-self.exprA.subexprs[k].left           
            my=self.exprB.subexprs[k].top -self.exprA.subexprs[k].top
            xa=mx*cubic(this_t)
            ya=my*cubic(this_t)           
            if(k=='term') or (k=='sign'):
                ya=ya-25
            if(k=='sign'):
                #Now we actually have to calculate the position of the other part
                fadea=1-this_t
                fadeb=this_t
                mx=self.exprB.subexprs[k].left-self.exprA.subexprs[k].left           
                my=self.exprB.subexprs[k].top -self.exprA.subexprs[k].top
                xb=mx*(cubic(this_t)-1)
                yb=my*(cubic(this_t)-1)-25           
        else:
            this_t=linterp(secondStageEnd,0,1,1,t)
            # No second part at all by default
            xb=0
            yb=0
            fadeb=0
            # first part completely visible
            fadea=1
            xa=self.exprB.subexprs[k].left-self.exprA.subexprs[k].left
            ya=self.exprB.subexprs[k].top -self.exprA.subexprs[k].top
            if(k=='term'):
                #move the term straight down
                ya+=cubic(1-this_t)*-25
            if(k=='sign'):
                #move the sign straight down and only show the second sign
                fadea=0
                fadeb=1
                xb=0
                yb=cubic(1-this_t)*-25
        return (xa,ya,fadea,xb,yb,fadeb)

class dance_transform_right(dance_eqn):
    def __init__(self,
                 left,
                 equals,
                 right_before,termA,termB,right_after):
        """
        Render the movement of a term from the left side to the right side.
    
        :param left:  Part of equation on the left side of the anchor in TeX format
        :type  left:  string
        :param equals:       Anchor of expression, usually equals. This will remain 
                             in the same place
        :type  equals:       string
        :param right_before: Part of equation on the right side to the left of 
                             (before) the term to be moved, in TeX format
        :type  right_before: string
        :param termA:   Term to be transformed, before it is transformed
        :type  termA:   string
        :param termB:   Term to be transformed, before it is transformed
        :type  termB:   string
        :param right_after:  Part of equation on the right side to the right of 
                             (after) the term to be moved
        :type  right_after:  string
        """
        # Get the spacing of the equation before the move
        self.exprA=expr("exprA",
                        ["left","equals","right_before","termA","right_after"],
                        [ left , equals , right_before , termA , right_after ],"equals")
        # Get the spacing of the equation after the move
        self.exprB=expr("exprB",
                        ["left","equals","right_before","termB","right_after"],
                        [ left , equals , right_before , termB , right_after ],"equals")
        
    def get_piece_choreography(self,k,t):
        firstStageEnd=0.2
        secondStageEnd=0.8
        xa=0
        ya=0
        fadea=1
        xb=0
        yb=0
        fadeb=1
        if(t<firstStageEnd):
            this_t=linterp(0,0,firstStageEnd,1,t)
            # No second part at all yet
            xb=0
            yb=0
            fadeb=0
            # first part completely visible
            fadea=1
            if(k=='term') or (k=='sign'):
                #move the first term and sign straight up
                xa=0
                ya=cubic(this_t)*-25
            else:
                xa=0
                ya=0
        elif(t<secondStageEnd):
            #No second part by default for most parts
            xb=0
            yb=0
            fadeb=0
            # first part completely visible
            fadea=1
            # Move parts by default
            this_t=linterp(firstStageEnd,0,secondStageEnd,1,t)
            mx=self.exprB.subexprs[k].left-self.exprA.subexprs[k].left           
            my=self.exprB.subexprs[k].top -self.exprA.subexprs[k].top
            xa=mx*cubic(this_t)
            ya=my*cubic(this_t)           
            if(k=='term') or (k=='sign'):
                ya=ya-25
            if(k=='sign'):
                #Now we actually have to calculate the position of the other part
                fadea=1-this_t
                fadeb=this_t
                mx=self.exprB.subexprs[k].left-self.exprA.subexprs[k].left           
                my=self.exprB.subexprs[k].top -self.exprA.subexprs[k].top
                xb=mx*(cubic(this_t)-1)
                yb=my*(cubic(this_t)-1)-25           
        else:
            this_t=linterp(secondStageEnd,0,1,1,t)
            # No second part at all by default
            xb=0
            yb=0
            fadeb=0
            # first part completely visible
            fadea=1
            xa=self.exprB.subexprs[k].left-self.exprA.subexprs[k].left
            ya=self.exprB.subexprs[k].top -self.exprA.subexprs[k].top
            if(k=='term'):
                #move the term straight down
                ya+=cubic(1-this_t)*-25
            if(k=='sign'):
                #move the sign straight down and only show the second sign
                fadea=0
                fadeb=1
                xb=0
                yb=cubic(1-this_t)*-25
        return (xa,ya,fadea,xb,yb,fadeb)

            
d1=dance_addsub_left_right("1","+","1","+\int_a^b x^2 dx","=","2","-","+\sum_{n=0}^\infty \sin(x^2) dx")
#d1=dance_transform_right(r"\frac{d}{dx}f(x)",r"=",r"\frac{",r"f(x+dx)",r"(x+dx)^2",r"-f(x)}{dx}")
fig,ax=plt.subplots(figsize=(16,9),dpi=108)
plt.axis('off')
frames=50
for i in range(frames):
    t=linterp(0,0,frames,1,i)
    print(t,i)
    plt.cla()
    ax.set_xlim(-200,200)
    ax.set_ylim(200*9/16,-200*9/16)
    d1.draw(ax,t,(0,0))
    plt.savefig("transform_%03d.png"%(i,))


