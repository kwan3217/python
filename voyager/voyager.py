import numpy as np
import matplotlib.pyplot as plt
import glob
import spiceypy as cspice


cspice.furnsh('/home/jeppesen/workspace/Data/spice/Voyager/sclk/vg100019.tsc')
cov0=cspice.ckcov('/home/jeppesen/workspace/Data/spice/Voyager/ck/vg1_jup_qmw_na.bc',-31001,False,"INTERVAL",0,"TDB")
print(cspice.wncard(cov0))

for infn in sorted(glob.glob('/home/jeppesen/workspace/Data/Voyager/Voyager 2 Io volcano Watch/*.IMG')):
    a=np.fromfile(infn,dtype='<i2')
    a.shape=(1001,1000)
    a=a[1:,:]
    print(infn,np.max(a))
    plt.figure(infn)
    plt.imshow(np.sqrt(np.maximum(a,0)/np.max(a)))
    plt.show()