import csv
import numpy as np
import matplotlib.pyplot as plt

pres_t=[]
pres=[]
mark_t=[]
mark=[]
at=[]
ax=[]
ay=[]
az=[]
gt=[]
gx=[]
gy=[]
gz=[]

with open('c:/Users/chrisj/Desktop/SensorLogs California 19.01/SpaceMtSensor.csv') as inf:
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
            mark.append(row[2])


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


mark_t=np.array(mark_t)
f_mark0=79
nf=12000
fdt=1001/30000 #Approximately 33ms, intended to be exact NTSC frame time ~29.97fps
ft=(np.arange(nf+1)-f_mark0)*fdt
print(ft[1],ft[f_mark0])
for this_mark_t,this_mark in zip(mark_t,mark):
    print((this_mark_t-mark_t[0])/1e9,this_mark)
pres_t=(np.array(pres_t,dtype=np.int64)-mark_t[0])/1e9
at=(np.array(at,dtype=np.int64)-mark_t[0])/1e9
ax=np.array(ax)
ay=np.array(ay)
az=np.array(az)
a=np.sqrt(ax**2+ay**2+az**2)
atrate=len(at)/(at[-1]-at[0])
print(atrate)
gt=(np.array(gt,dtype=np.int64)-mark_t[0])/1e9
gx=np.degrees(np.array(gx))
gy=np.degrees(np.array(gy))
gz=np.degrees(np.array(gz))
gtrate=len(gt)/(gt[-1]-gt[0])
print(gtrate)
pres=np.array(pres)
pres_alt=145366.45*(1-(pres/1013.25)**0.190284)
N=int(atrate*0.25)
mat=running_mean(at,N)
ma=running_mean(a,N)
mpres_t=running_mean(pres_t,N)
mpres_alt=running_mean(pres_alt,N)
fig=plt.figure(figsize=[16,9])
ax=fig.add_subplot(111)
pag,=ax.plot(at, a, color="#c0c0c0")
ppg,=ax.plot(pres_t, pres_alt, color="#c0c0c0")
pa,=ax.plot(mat, ma, 'r-')
pp,=ax.plot(mpres_t, running_mean(pres_alt, N), 'k-')
ph0,=ax.plot([0,0],[-50,50],'k-')
phm,=ax.plot([0,0],[-50,50],color="#808080")
php,=ax.plot([0,0],[-50,50],color="#808080")
for i in range(1,nf):
    pag.set_xdata(at)
    pag.set_ydata(a)
    ppg.set_xdata(pres_t)
    ppg.set_ydata(pres_alt)
    pa.set_xdata(mat)
    pa.set_ydata(ma)
    pp.set_xdata(mpres_t)
    pp.set_ydata(mpres_alt)
    ph0.set_xdata([ft[i],ft[i]])
    phm.set_xdata([ft[i]-0.125,ft[i]-0.125])
    php.set_xdata([ft[i]+0.125,ft[i]+0.125])
    plt.ylim(-50,50)
    plt.xlim(ft[i]-7,ft[i]+3)
    fig.savefig("Frames/frame%05d.png"%i,dpi=120,transparent=True)
    plt.pause(0.01)