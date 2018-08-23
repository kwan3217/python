import glob
import struct
import os.path
import numpy as np
import matplotlib.pyplot as plt

chunks=0
regions=0
chunkdict = {}
regiondict={}
cxmax = float('-inf')
cxmin = float('inf')
czmax = float('-inf')
czmin = float('inf')
rxmax = float('-inf')
rxmin = float('inf')
rzmax = float('-inf')
rzmin = float('inf')
for infn in glob.glob("/home/jeppesen/Desktop/backup20180823/region/r.*.*.mca"):
    base_infn=os.path.basename(infn)
    (_,rx,rz,_)=base_infn.split(".")
    rx=int(rx)
    rz=int(rz)
    regiondict[(rx,rz)]=True
    regions+=1
    if rx > rxmax:
        rxmax = rx
    if rx < rxmin:
        rxmin = rx
    if rz > rzmax:
        rzmax = rz
    if rz < rzmin:
        rzmin = rz
    print(infn)
    with open(infn,"rb") as inf:
        header=inf.read(4096)
        offsets=struct.unpack(">"+"L"*1024,header)
        for i,offset in enumerate(offsets):
            cx=rx*32+i % 32
            cz=rz*32+i // 32
            if offset>0:
                if cx>cxmax:
                    cxmax=cx
                elif cx<cxmin:
                    cxmin=cx
                if cz>czmax:
                    czmax=cz
                elif cz<czmin:
                    czmin=cz
                chunkdict[(cx,cz)]=True
                chunks+=1

def dict_to_map(dict,xmin,xmax,zmin,zmax):
    map=np.zeros((zmax-zmin+1,xmax-xmin+1))
    for x,z in dict:
        map[z-zmin,x-xmin]=1
    map[-zmin,-xmin]=2
    map[-zmin-1,-xmin]=3 #arrow pointing north
    return map

regionmap=dict_to_map(regiondict,rxmin,rxmax,rzmin,rzmax)
print(chunks,cxmin,cxmax,czmin,czmax)
plt.imshow(regionmap,origin='upper',extent=(rxmin,rxmax+1,rzmax+1,rzmin))
plt.show()
chunkmap=dict_to_map(chunkdict,cxmin,cxmax,czmin,czmax)
plt.imshow(chunkmap,origin='upper',extent=(cxmin*16,cxmax*16+16,czmax*16+16,czmin*16))
plt.show()