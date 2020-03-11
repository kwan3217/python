import glob
import struct
import os.path
import numpy as np
import matplotlib.pyplot as plt
import nbt
import zlib
import io

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
for infn in glob.glob("/home/minecraft/world-1.15.2-20200310/region/r.*.*.mca"):
    base_infn=os.path.basename(infn)
    print(base_infn)
    (_,rx,rz,_)=base_infn.split(".")
    rx=int(rx)
    rz=int(rz)
    print(infn)
    with open(infn,"rb") as inf:
        header=inf.read(4096)
        offsets=struct.unpack(">"+"L"*1024,header)
        region_chunks=0
        for i,offset in enumerate(offsets):
            cx=rx*32+i % 32
            cz=rz*32+i // 32
            if offset>0:
                sector_ofs = (offset >> 8 & 0x00FFFFFF)*0x1000
                len_sectors = (offset & 0xFF)*0x1000
                inf.seek(sector_ofs)
                buf=inf.read(5)
                (exact_len,compression)=struct.unpack(">LB",buf)
                buf=inf.read(exact_len)
                buf=zlib.decompress(buf,47)
                chunk=nbt.nbt.NBTFile(buffer=io.BytesIO(buf))
                region_chunks += 1
                if cx>cxmax:
                    cxmax=cx
                elif cx<cxmin:
                    cxmin=cx
                if cz>czmax:
                    czmax=cz
                elif cz<czmin:
                    czmin=cz
                chunkdict[(cx,cz)]=chunk['Level']['InhabitedTime'].value
                chunks+=1
    if region_chunks>0:
        regiondict[(rx,rz)]=1
        regions+=1
        if rx > rxmax:
            rxmax = rx
        if rx < rxmin:
            rxmin = rx
        if rz > rzmax:
            rzmax = rz
        if rz < rzmin:
            rzmin = rz

def dict_to_map(dict,xmin,xmax,zmin,zmax,mark=True):
    map=np.zeros((zmax-zmin+1,xmax-xmin+1))
    for ((x,z),value) in dict.items():
        map[z-zmin,x-xmin]=value
    map[-zmin,-xmin]=2
    map[-zmin-1,-xmin]=3 #arrow pointing north
    return map

regionmap=dict_to_map(regiondict,rxmin,rxmax,rzmin,rzmax)
print("Regions: ",regions,rxmin,rxmax,rzmin,rzmax)
plt.imshow(regionmap,origin='upper',extent=(rxmin,rxmax+1,rzmax+1,rzmin))
plt.show()
chunkmap=dict_to_map(chunkdict,cxmin,cxmax,czmin,czmax)
print("Chunks: ",chunks,cxmin,cxmax,czmin,czmax)
plt.imshow(np.log2(chunkmap+1),origin='upper',extent=(cxmin*16,cxmax*16+16,czmax*16+16,czmin*16))
plt.show()

