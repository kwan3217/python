import csv
import numpy as np
import matplotlib.pyplot as plt

pres_t=[]
pres=[]
mark_t=[]
at=[]
ax=[]
ay=[]
az=[]
gt=[]
gx=[]
gy=[]
gz=[]

with open('/home/chrisj/Desktop/SensorLogs California 19.01/SpaceMtSensor.csv') as inf:
    csvreader=csv.reader(inf,delimiter=',')
    for row in csvreader:
        if row[1]=="1":
            #Calibrated accelerometer measurement
            at.append(int(row[0]))
            ax.append(float(row[2]))
            ay.append(float(row[3]))
            az.append(float(row[4]))
        if row[1]=="4":
            #Calibrated gyroscope measurement
            gt.append(int(row[0]))
            gx.append(float(row[2]))
            gy.append(float(row[3]))
            gz.append(float(row[4]))
        if row[1]=="6":
            #Pressure measurement
            pres_t.append(int(row[0]))
            pres.append(float(row[2]))
        if row[1]=="65537":
            mark_t.append(int(row[0]))


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


pres_t=(np.array(pres_t,dtype=np.int64)-mark_t[0])/1e9
at=(np.array(at,dtype=np.int64)-mark_t[0])/1e9
ax=np.array(ax)
ay=np.array(ay)
az=np.array(az)
a=np.sqrt(ax**2+ay**2+az**2)
gt=(np.array(gt,dtype=np.int64)-mark_t[0])/1e9
gx=np.degrees(np.array(gx))
gy=np.degrees(np.array(gy))
gz=np.degrees(np.array(gz))
pres=np.array(pres)
pres_alt=145366.45*(1-(pres/1013.25)**0.190284)
N=500
plt.plot(running_mean(at,N),running_mean(a,N),'r-',running_mean(pres_t,N),running_mean(pres_alt,N),'k-')
plt.show()
