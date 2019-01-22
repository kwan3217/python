'''
Created on Aug 3, 2017

@author: chrisj
'''
import glob
import struct
import xtce_packet

def parseFile(infn,pktDefs):
    basefn=".".join(infn.split(".")[0:-1])
    inf=open(infn,"rb")
    inf.read(8)
    tmin=[0]*64
    tminlast=tmin[:]
    tlastfast=[0]*64
    ouf=[None]*64
    oufName=[None,                #Apid 0
             None,                #Apid 1
             (basefn+".bin"      ,"wb"), #Apid 2
             (basefn+".cpio.zpaq","wb"), #Apid 3
             None,                #Apid 4
             None,                #Apid 5
             None,                #Apid 6
             None,                #Apid 7
             None,                #Apid 8
             (basefn+".%02d%02d.mpu.csv","w"), #Apid  9
             (basefn+".%02d%02d.mag.csv","w"), #Apid 10
             (basefn+".bmp.csv","w"), #Apid 11
             None,                #Apid 12
             (basefn+".nmea"            ,"w"), #Apid 13
             None,                #Apid 14
             None]                #Apid 15
    pclk=0
    old_temp_cal=0
    old_pres_cal=0
    while True:
        try:
            packet=xtce_packet.getNextPacket(inf,pktDefs,[1,11,13])
        except struct.error:
            break
            for f in ouf:
                f.close()
        apid=packet["PKT_APID"]
        if "TC" in packet:
            if packet["TC"]<tlastfast[apid]:
                if oufName[apid][0].find("%")>=0 and ouf[apid] is not None:
                    ouf[apid].close()
                    ouf[apid]=None
                tmin[apid]+=1
            tlastfast[apid]=packet["TC"]
            packet["TC"]="%02d:%02d:%011.8f"%(tmin[apid]//60,tmin[apid]%60,float(packet["TC"])/pclk)
        if ouf[apid] is None:
            if oufName[apid] is not None:
                if oufName[apid][0].find("%")>=0:
                    ouf[apid]=open(oufName[apid][0]%(tmin[apid]//60,tmin[apid]%60),oufName[apid][1])
                else:
                    ouf[apid]=open(oufName[apid][0],oufName[apid][1])
        if apid==1:
            pclk=float(packet["PCLK"])
        elif apid in [2,3]:
            ouf[apid].write(packet["DATA"])
        elif packet["PKT_APID"]==9:
            ouf[apid].write("%s,%d,%d,%d,%d,%d,%d,%d\n"%(packet["TC"],packet["AX"],packet["AY"],packet["AZ"],packet["GX"],packet["GY"],packet["GZ"],packet["MT"]))
        elif packet["PKT_APID"]==10:
            ouf[apid].write("%s,%d,%d,%d\n"%(packet["TC"],packet["BX"],packet["BY"],packet["BZ"]))
        elif packet["PKT_APID"]==11:
            if packet["TEMP_CAL"]!=old_temp_cal or packet["PRES_CAL"]!=old_pres_cal:
                ouf[apid].write("%s,%4.1f,%d\n"%(packet["TC"],packet["TEMP_CAL"],packet["PRES_CAL"]))
                old_temp_cal=packet["TEMP_CAL"]
                old_pres_cal=packet["PRES_CAL"]
        elif packet["PKT_APID"]==13:
            ouf[apid].write(packet["NMEA"])
        else:
            print(packet)

if __name__ == '__main__':
    pktDefs=xtce_packet.readPktDesc('xtce_pocketometer20150317.xml')
    for infn in glob.glob("*.SDS"):
        print(infn)
        parseFile(infn,pktDefs)
