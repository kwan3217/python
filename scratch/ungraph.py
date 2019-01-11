import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import scipy.interpolate

#Screen coordinates of axes in columns. 1 at bottom row of matrix is required to make it square
X=np.array([[326.75,1913.5,326.5],
            [990,989.5,81.5],
            [1,1,1]])
print("X: ",X)
#Data coordinates of latitude axes in columns
Blat=np.array([[34680,35520,34680],
            [   24,   24,   36],
            [    1,    1,    1]])
print("Blat: ",Blat)
#Data coordinates of east longitude axes in columns
Blon=np.array([[34680,35520,34680],
               [  156,  156,  204],
               [    1,    1,    1]])
print("Blon: ",Blon)
#Data coordinates of altitude axes in columns
Balt=np.array([[34680,35520,34680],
               [-40e3,-40e3,440e3],
               [    1,    1,    1]])
print("Balt: ",Balt)
#Transformation from pixel to data coordinates
#A@X=B
#A@X@X^-1=B@X^-1
#A       =B@X^-1
Alat=Blat @ numpy.linalg.inv(X)
print("Alat: ",Alat)
Alon=Blon @ numpy.linalg.inv(X)
print("Alon: ",Alon)
Aalt=Balt @ numpy.linalg.inv(X)
print("Aalt: ",Aalt)
#Test - all elements should be very small, near machine epsilon, ideally zero
print("lat test: ",Blat-(Alat @ X))
print("lon test: ",Blon-(Alon @ X))
print("alt test: ",Balt-(Aalt @ X))

#Time scale
TimeData=[[364,989.75],
[401.5,990.25],
[439.5,990.75],
[476.5,990.5],
[515.75,990.5],
[553.25,990.5],
[591,990.5],
[628.25,990.75],
[665.75,990.5],
[703.5,990.5],
[741.5,990.5],
[779,990.5],
[817,990.5],
[854.5,990.5],
[891.5,990.75],
[929.5,990.75],
[967,990.25],
[1004,990.5],
[1042.25,990.75],
[1080.5,990.75],
[1118.75,990.5],
[1156.75,990.25],
[1194.5,990.5],
[1232.25,990.5],
[1270,990.5],
[1307.75,990.5],
[1345.25,990.5],
[1383,990],
[1421.5,990],
[1458.5,990],
[1496.5,989.75],
[1534,990],
[1572.25,990.25],
[1609.5,990],
[1648,989.75],
[1685,989.75],
[1723.5,989.5],
[1761.5,989.25],
[1799.5,989.5],
[1838,989.25],
[1875.5,989.5]]

TimeScale=np.zeros((3,len(TimeData)+2))
TimeScale[:,0]=X[:,0]
for i in range(len(TimeData)):
    TimeScale[:,i+1]=np.array([TimeData[i]+[1]])
TimeScale[:,-1]=X[:,1]
print("TimeScale: ",TimeScale)
print("TimeDataScale: ",Alat @ TimeScale)
TimeData=(Alat @ TimeScale)[0,:]
print("TimeData: ",TimeData)
RoundedTimeData=np.round(TimeData/20)*20
print("RoundedTimeData: ",RoundedTimeData)
TimeMap={}
ReverseTimeMap={}
for i in range(len(TimeData)):
    TimeMap[TimeData[i]]=RoundedTimeData[i]
    ReverseTimeMap[RoundedTimeData[i]]=TimeData[i]
print("TimeMap: ",TimeMap)
print("ReverseTimeMap: ",ReverseTimeMap)
#plt.plot(RoundedTimeData,RoundedTimeData-TimeData)
#plt.show()

#Vertical grid scale
VertData=[[326.75,952.25],
[327,914.5],
[327,876],
[327,838.5],
[327.25,801],
[327,763],
[327,725.25],
[327.25,686.75],
[326.75,649.25],
[326.75,611.75],
[327,574.25],
[326.75,536],
[327,498.5],
[326.75,460.25],
[326.75,422],
[326.5,384.75],
[326.5,346.5],
[326.5,308.75],
[326.5,271.5],
[326.5,233],
[326.5,195.5],
[326.5,157.5],
[326.5,119.25]]

VertScale=np.zeros((3,len(VertData)+2))
VertScale[:,0]=X[:,0]
for i in range(len(VertData)):
    VertScale[:,i+1]=np.array([VertData[i]+[1]])
VertScale[:,-1]=X[:,2]
print("VertScale: ",VertScale)
print("Alat @ VertScale: ",Alat @ VertScale)
VertData=(Alat @ VertScale)[1,:]
print("VertData: ",VertData)
RoundedVertData=np.round(VertData/0.5)*0.5
print("RoundedVertData: ",RoundedVertData)
#plt.plot(RoundedVertData,RoundedVertData-VertData)
#plt.show()

#Latitude data
LatData=Alat @ np.array([[379,329,1],
[401,330,1],
[439.5,333,1],
[476.75,337.75,1],
[515.25,345.25,1],
[552.25,354.25,1],
[590.5,363.5,1],
[627.75,372.5,1],
[665.25,381,1],
[703.25,390.25,1],
[741,400.75,1],
[778.5,412,1],
[816.25,426,1],
[853.75,441,1],
[891.25,457,1],
[929,473.5,1],
[966.25,491,1],
[1002.75,509.5,1],
[1041.5,527.75,1],
[1079.25,546.75,1],
[1117.5,566.75,1],
[1156,585.5,1],
[1193.5,602.5,1],
[1231.75,618.75,1],
[1269,634.75,1],
[1307.25,647.25,1],
[1344.5,659.75,1],
[1382.5,672,1],
[1420.5,682,1],
[1457.75,692.5,1],
[1496.25,698.75,1],
[1533.75,704.5,1],
[1571.5,709,1],
[1609,712,1],
[1647,715,1],
[1684.75,716.75,1],
[1722.75,718,1],
[1760.75,718.5,1],
[1798.75,719.25,1],
[1837,719.5,1],
[1850.25,719.75,1]]).T

def FixTimeSeries(Data):
    global TimeData,TimeMap,ReverseTimeMap
    for i in range(Data.shape[1]):
        rounded=np.round(Data[0,i]/20)*20
        if np.abs(Data[0,i]-ReverseTimeMap[rounded])<3:
            #Grid snap
            Data[0,i]=rounded
        else:
            #Figure out the linear interpolated time mapping
            #which segment are we in?
            for j in range(len(TimeData)):
                if Data[0,i]<TimeData[j]:
                    t=(Data[0,i]-TimeData[j-1])/(TimeData[j]-TimeData[j-1])
                    Data[0,i]=(1-t)*TimeMap[TimeData[j-1]]+t*TimeMap[TimeData[j]]
                    break

FixTimeSeries(LatData)
print("LatData: ",LatData)

#Latitude data
LonData=Alon @ np.array([[378.75,799,1],
[401.5,776.75,1],
[439.75,736.75,1],
[476.75,700.25,1],
[516,660.75,1],
[553,627,1],
[590.75,593.25,1],
[628.75,561.75,1],
[666,532.25,1],
[704,504.75,1],
[741,477.75,1],
[778.5,453.25,1],
[816.25,429,1],
[853.5,405.5,1],
[891.5,381.25,1],
[928.25,359.5,1],
[966,336.75,1],
[1002,315,1],
[1041.5,294.75,1],
[1078.75,274.75,1],
[1117,257,1],
[1154.75,239.25,1],
[1193,222.25,1],
[1231.25,207,1],
[1268.25,194.25,1],
[1306,181.75,1],
[1342.75,171.5,1],
[1381.75,161.25,1],
[1418.75,151.5,1],
[1456.75,143,1],
[1494,137.5,1],
[1532.5,133,1],
[1570,129,1],
[1607.5,125.75,1],
[1645.25,124,1],
[1682.75,122.75,1],
[1721,122,1],
[1758.25,121.5,1],
[1796.5,121.5,1],
[1834,121.25,1],
[1848.25,121,1]]).T
FixTimeSeries(LonData)
print("LonData: ",LonData)
#plt.plot(LonData[1,:],LatData[1,:],'*-b')
#plt.show()

#Latitude data
AltData=Aalt @ np.array([[375.75,158,1],
[401.25,229.5,1],
[439.5,336.25,1],
[477,431.25,1],
[516,505.5,1],
[535.75,536,1],
[552.25,552,1],
[560,556.25,1],
[568.25,558.25,1],
[576.75,557,1],
[585,555,1],
[591.25,552.75,1],
[610,545.25,1],
[628.75,538.75,1],
[665.75,528.25,1],
[703.75,522.75,1],
[741.5,518.75,1],
[778.5,513.5,1],
[816.75,507.25,1],
[854,502.5,1],
[872.5,501.5,1],
[891.25,502.5,1],
[928.75,507.25,1],
[966.5,514.75,1],
[1003,525.25,1],
[1042.25,541,1],
[1079.25,556.75,1],
[1117.5,571,1],
[1155.75,582.75,1],
[1193.25,591.75,1],
[1231.75,600.25,1],
[1269.25,609.5,1],
[1307.25,619.25,1],
[1344.5,630,1],
[1382.5,642.5,1],
[1420.75,657.25,1],
[1457.75,672.75,1],
[1496.25,688.75,1],
[1534.25,706,1],
[1571.5,723.75,1],
[1609.25,742.5,1],
[1647.5,763,1],
[1684.5,783.75,1],
[1723.25,806.5,1],
[1761,827.5,1],
[1799,849.25,1],
[1829.5,866.5,1]]).T
FixTimeSeries(AltData)
print("AltData: ",AltData)
plt.plot(AltData[0,:],AltData[1,:],'*-b')
plt.show()

t0=np.min([LatData[0,0],LonData[0,0],AltData[0,0]])
t1=np.max([LatData[0,-1],LonData[0,-1],AltData[0,-1]])
LatData[0,0]=t0
LonData[0,0]=t0
AltData[0,0]=t0
LatData[0,-1]=t1
LonData[0,-1]=t1
AltData[0,-1]=t1

Lat=scipy.interpolate.interp1d(LatData[0,:],LatData[1,:],kind='cubic')
Lon=scipy.interpolate.interp1d(LonData[0,:],LonData[1,:],kind='cubic')
Alt=scipy.interpolate.interp1d(AltData[0,:],AltData[1,:],kind='cubic')

print("t,lon (degE),lat (deg),alt (ft)")
for t in AltData[0,:]:
    print(t,",",Lon(t),",",Lat(t),",",Alt(t))

plt.plot(Lon(AltData[0,:]),Alt(AltData[0,:]),"*-b")
plt.show()