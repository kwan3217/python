import glob
import sys
import os.path

for infn in glob.glob("UltimateGPS/*/*.txt"):
    oufn=os.path.basename(infn)
    oufn=oufn.split(".")
    oufn=".".join(oufn[:-1])
    oufn=os.path.dirname(infn)+"/"+oufn+".nmea"
    print(infn)
    print(oufn)
    with open(infn,"r") as inf:
        with open(oufn,"w") as ouf:
            for line in inf:
                if line[0]=="$":
                    ouf.write(line)

