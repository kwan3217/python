import subprocess
import os

def tex_head(oufn):
    ouf=open(oufn+'.tex','wt')
    print(r"\documentclass{minimal} \batchmode",file=ouf)
    print(r"\usepackage[usenames,dvipsnames,svgnames,table]{xcolor}",file=ouf)
    print(r"\begin{document}",file=ouf)
    print(r"$$",file=ouf)
    return ouf

def tex_eqn(ouf,text):
    print("{"+text+"}",file=ouf)

def tex_foot(ouf):
    print("$$",file=ouf)
    print("\end{document}",file=ouf)
    ouf.close()

def tex_render(oufn,clean=False,verbose=True):
    # Get the coordinates of the piece by itself
    subprocess.call("latex "+oufn+".tex"+(" > /dev/null" if not verbose else ""),shell=True)
    if clean:
        # Delete the .tex since we no longer need it
        os.remove(oufn+'.tex')
        # clean up other files we don't care about
        os.remove(oufn+'.aux')
        os.remove(oufn+'.log')
    # Render the piece
    subprocess.call("dvisvgm -e -n -bmin --keep "+oufn+".dvi > /dev/null 2>&1",shell=True)
    # Delete the DVI now that it is rendered
    if clean:
        os.remove(oufn+'.dvi')
    # Slurp the SVG for this equation
    with open(oufn+'.svg',"r") as inf:
        content=inf.readlines()
    # Delete the SVG
    if clean:
        os.remove(oufn+'.svg')
    return content

def eqn2svg(eqn,eqnname="temp"):
    ouf=tex_head(eqnname)
    tex_eqn(ouf,eqn)
    tex_foot(ouf)
    return tex_render(eqnname)

def test_eqn2svg():
    print(eqn2svg(r"y=\frac{-c \pm \sqrt{d^2-{\color[rgb]{1,0,0}4ef}}}{2g}"))

if __name__=="__main__":
    test_eqn2svg()