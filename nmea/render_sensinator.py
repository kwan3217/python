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

with open('SpaceMtSensor.csv') as inf:
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
gx=np.array(gx)
gy=np.array(gy)
gz=np.array(gz)
g=np.sqrt(gx**2+gy**2+gz**2)
gtrate=len(gt)/(gt[-1]-gt[0])
print(gtrate)
pres=np.array(pres)
pres_alt=145366.45*(1-(pres/1013.25)**0.190284)
pres_trate=len(pres_t)/(pres_t[-1]-pres_t[0])
N=int(atrate*0.25)
mat=running_mean(at,N)
ma=running_mean(a,N)
N=int(gtrate*0.25)
mgt=running_mean(gt,N)
mg=running_mean(g,N)
N=int(pres_trate*0.25)
mpres_t=running_mean(pres_t,N)
mpres_alt=running_mean(pres_alt,N)
fig=plt.figure(figsize=[16,9])
ax=fig.add_subplot(111)
ph0,=ax.plot([0,0],[-50,50],color="#ffffff")
phm,=ax.plot([0,0],[-50,50],color="#e0e0e0")
php,=ax.plot([0,0],[-50,50],color="#e0e0e0")
pam,=ax.plot(at, a, color="#800000")
pgm,=ax.plot(gt, g, color="#008000")
ppm,=ax.plot(pres_t, pres_alt, color="#000080")
pap,=ax.plot(at, a, color="#808080")
pgp,=ax.plot(gt, g, color="#808080")
ppp,=ax.plot(pres_t, pres_alt, color="#808080")
pa,=ax.plot(mat, ma, color="#ff0000")
pg,=ax.plot(mgt, mg, color="#00ff00")
pp,=ax.plot(mpres_t, mpres_alt,color="#0000ff")
an=ax.text(0,0,"acc: %4.2g"%0,color="#ff0000")
gn=ax.text(0,0,"rot: %4.2rad/s"%0,color="#00ff00")
pn=ax.text(0,0,"palt: %4.1ft"%0,color="#0000ff")
#https://stackoverflow.com/a/12059429
#When using figures, you can easily change the spine color with:
ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')
ax.spines['left'].set_color('white')
#Use the following to change only the ticks:
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
#And the following to change only the label:
ax.yaxis.label.set_color('white')
ax.xaxis.label.set_color('white')
#And finally the title:
ax.title.set_color('white')

def update_line(traw,raw,tsmooth,smooth,ft,line_smooth,line_rawm,line_rawp,text,label,scale=1.0):
    i=np.abs(traw-ft).argmin()
    ismooth=np.abs(tsmooth-ft).argmin()
    text.set_text(label%(smooth[ismooth]/scale))
    text.set_x(ft)
    text.set_y(smooth[ismooth]+1)
    line_rawm.set_xdata(traw[:i])
    line_rawm.set_ydata(raw[:i])
    line_rawp.set_xdata(traw[i:])
    line_rawp.set_ydata(raw[i:])
    line_smooth.set_xdata(tsmooth[:ismooth])
    line_smooth.set_ydata(smooth[:ismooth])

for i in range(1,nf):
    plt.ylim(-50,50)
    plt.xlim(ft[i]-7,ft[i]+3)
    update_line(at,a,mat,ma,ft[i],pa,pam,pap,an,"acc: %4.2fg",9.80665)
    update_line(gt,g,mgt,mg,ft[i],pg,pgm,pgp,gn,"rot: %4.2frad/s")
    update_line(pres_t,pres_alt,mpres_t,mpres_alt,ft[i],pp,ppm,ppp,pn,"alt: %4.2fft")
    print(i)
    if i==289:
        print(i)
    ph0.set_xdata([ft[i],ft[i]])
    phm.set_xdata([ft[i]-0.125,ft[i]-0.125])
    php.set_xdata([ft[i]+0.125,ft[i]+0.125])
    plt.pause(0.01)
    fig.savefig("Frames/frame%05d.png"%i,dpi=120,transparent=True)
