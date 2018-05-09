"""
Code for parsing Spice Double-Array Files (DAF). This is a superset of CK and SPK files. Note that while this
parses Spice files, it only makes minimal use of the spice libraries (date conversion only)
"""

import struct
from spiceypy import etcal  as cspice_etcal
from spiceypy import scdecd as cspice_scdecd
from spiceypy import sct2e  as cspice_sct2e

class daf_summary:
    @staticmethod
    def make(sr,daf,buf,name):
        if daf.subtype=="SPK":
            return daf_SPKsummary(sr,daf,buf,name)
        elif daf.subtype=="CK":
            return daf_CKsummary(sr,daf,buf,name)
        else:
            raise ValueError("Unrecognized subtype %s"%daf.subtype)
    def __init__(self,sr,daf,buf,name):
        """

        :param sr: Summary record which this summary is part of
        :param daf: Double-array file which holds the header (ND, NI) and can read the file
        :param buf: Data for this summary
        """
        self.sr=sr
        self.daf=daf
        self.D=[0.0]*self.daf.ND
        self.I=[0]*self.daf.NI
        self.name=name
        for i in range(self.daf.ND):
            self.D[i]=struct.unpack(self.daf.dformat,buf[i*8:(i+1)*8])[0]
        for i in range(self.daf.NI):
            self.I[i]=struct.unpack(self.daf.iformat,buf[self.daf.ND*8+i*4:self.daf.ND*8+(i+1)*4])[0]
    def __str__(self):
        result="%s\nD:\n" % self.name
        for i,D in enumerate(self.D):
            if i!=0:
                result+=","
            result+="%20.15e"%D
        result+="\nI:\n"
        for i,I in enumerate(self.I):
            if i!=0:
                result+=","
            result+="%d"%I
        return result

class daf_SPKsummary(daf_summary):
    def __init__(self,sr,daf,buf,name):
        super().__init__(sr,daf,buf,name)
    def getet0(self): return self.D[0]
    def getet1(self): return self.D[1]
    def gettarget(self): return self.I[0]
    def getcenter(self): return self.I[1]
    def getframe (self): return self.I[2]
    def gettype  (self): return self.I[3]
    def getaddr0 (self): return self.I[4]
    def getaddr1 (self): return self.I[5]
    et0   =property(fget=getet0   ,doc="Ephemeris time of beginning of segment")
    et1   =property(fget=getet1   ,doc="Ephemeris time of end of segment")
    target=property(fget=gettarget,doc="Target NAIF code")
    center=property(fget=getcenter,doc="Center NAIF code")
    frame =property(fget=getframe ,doc="Reference frame NAIF code")
    type  =property(fget=gettype  ,doc="SPK segment type")
    addr0 =property(fget=getaddr0 ,doc="Address of first element of segment data")
    addr1 =property(fget=getaddr1 ,doc="Address of last element of segment data")
    def __str__(self):
        result=(("%s\n"+
                 "ET0:    %30.14f (%s)\n"+
                 "ET1:    %30.14f (%s)\n"+
                 "Target: %d\n"+
                 "Center: %d\n"+
                 "Frame:  %d\n"+
                 "Type:   %d\n"+
                 "Addr0:  %d\n"+
                 "Addr1:  %d") %
                (self.name,self.et0,cspice_etcal(self.et0),self.et1,cspice_etcal(self.et1),self.target,self.center,self.frame,self.type,self.addr0,self.addr1))
        return result
    def segment(self):
        return daf_SPKSegment(self,self.daf)

class daf_CKsummary(daf_summary):
    def __init__(self,sr,daf,buf,name):
        super().__init__(sr,daf,buf,name)
    def getsclk0(self): return self.D[0]
    def getsclk1(self): return self.D[1]
    def gettarget(self): return self.I[0]
    def getframe (self): return self.I[1]
    def gettype  (self): return self.I[2]
    def getrates (self): return self.I[3]
    def getaddr0 (self): return self.I[4]
    def getaddr1 (self): return self.I[5]
    sclk0 =property(fget=getsclk0 ,doc="DP encoded sclk of beginning of segment")
    sclk1 =property(fget=getsclk1 ,doc="DP encoded sclk of end of segment")
    target=property(fget=gettarget,doc="Target NAIF code")
    frame =property(fget=getframe ,doc="Reference frame NAIF code")
    type  =property(fget=gettype  ,doc="CK segment type")
    rates =property(fget=getrates ,doc="Angular Rates Flag")
    addr0 =property(fget=getaddr0 ,doc="Address of first element of segment data")
    addr1 =property(fget=getaddr1 ,doc="Address of last element of segment data")
    def __str__(self):
        result=(("%s\n"+
                 "SCLK0:    %30.14f\n"+
                 "SCLK1:    %30.14f\n"+
                 "Instrument: %d\n"+
                 "Frame:  %d\n"+
                 "Type:   %d\n"+
                 "Rates?  %d\n"+
                 "Addr0:  %d\n"+
                 "Addr1:  %d") %
                (self.name,self.sclk0,self.sclk1,self.target,self.frame,self.type,self.rates,self.addr0,self.addr1))
        return result
    def segment(self):
        return daf_CKSegment(self,self.daf)

class daf_SPKSegment:
    csvheaders=[0]*14
    csvheaders[1]="ET,ETCAL*"
    for i in range(71):
        csvheaders[1]+=",e_%02d"%i
    csvheaders[9]="ET,ETCAL*,x,y,z,dxdt,dydt,dzdt"
    csvheaders[13]=csvheaders[9]
    def __init__(self,summary,daf):
        self.summary=summary
        self.daf=daf
        daf.inf.seek(8*(summary.addr0-1))
        segment_size=summary.addr1-summary.addr0+1
        format="%s%d%s"%(daf.dformat[0],segment_size,daf.dformat[1])
        self.buf=struct.unpack(format,daf.inf.read(segment_size*8))
        self.N=int(self.buf[-1])
    def line(self,i):
        if self.summary.type==1:
            return daf_SPK01line(self.N,i,self.buf)
        elif self.summary.type==9 or self.summary.type==13:
            return daf_SPK09line(self.N,i,self.buf)
    def __iter__(self):
        for i in range(self.N):
            result = self.line(i)
            yield result
    def csvheader(self):
        return self.csvheaders[self.summary.type]
    def __str__(self):
        result=(("N: %d") %
                (self.N))
        return result

class daf_CKSegment:
    csvheaders=[0]*14
    csvheaders[3]="SCLK,q0,q1,q2,q3,wx,wy,wz"
    def __init__(self,summary,daf):
        self.summary=summary
        self.daf=daf
        daf.inf.seek(8*(summary.addr0-1))
        segment_size=summary.addr1-summary.addr0+1
        format="%s%d%s"%(daf.dformat[0],segment_size,daf.dformat[1])
        self.buf=struct.unpack(format,daf.inf.read(segment_size*8))
        self.N=int(self.buf[-1])
    def line(self,i):
        if self.summary.type==3:
            return daf_CK03line(self.N,i,self.buf,self.summary.rates)
        elif self.summary.type==9 or self.summary.type==13:
            return daf_SPK09line(self.N,i,self.buf)
    def __iter__(self):
        for i in range(self.N):
            result = self.line(i)
            yield result
    def csvheader(self):
        return self.csvheaders[self.summary.type]
    def __str__(self):
        result=(("N: %d") %
                (self.N))
        return result


class daf_SPK09line:
    def __init__(self,N,i,buf):
        self.x   =buf[i*6+0]
        self.y   =buf[i*6+1]
        self.z   =buf[i*6+2]
        self.dxdt=buf[i*6+3]
        self.dydt=buf[i*6+4]
        self.dzdt=buf[i*6+5]
        self.et=buf[N*6+i]
    def __str__(self):
        result="%30.14f,%s,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e"%(self.et,cspice_etcal(self.et),self.x,self.y,self.z,self.dxdt,self.dydt,self.dzdt)
        return result

class daf_SPK01line:
    n_elements=71
    def __init__(self,N,i,buf):
        self.e   =buf[i*self.n_elements:(i+1)*self.n_elements]
        self.et=buf[N*self.n_elements+i]
    def __str__(self):
        result="%30.14f,%s" %(self.et,cspice_etcal(self.et))
        for i in range(self.n_elements):
            result+=",%25.15e"%self.e[i]
        return result

class daf_CK03line:
    def __init__(self,N,i,buf,rates):
        self.rates=rates
        if self.rates:
            self.n_elements = 7
            self.w = buf[i * self.n_elements+4:(i + 1) * self.n_elements]
        else:
            self.n_elements=4
        self.q = buf[i * self.n_elements:i * self.n_elements+4]

        self.sclk=buf[N*self.n_elements+i]
    def __str__(self):
        et=cspice_sct2e(-84, self.sclk)
        result="%30.14f,%s,%30.14f,%s" %(self.sclk,cspice_scdecd(-84,self.sclk),et,cspice_etcal(et))
        for i in range(4):
            result+=",%25.15e"%self.q[i]
        if self.rates:
            for i in range(3):
                result += ",%25.15e" % self.w[i]
        return result

class daf_summary_record:
    def __init__(self,parent):
        """
        Read a record from a file-like object opened in binary mode
        :param parent: object of class daf
        """
        self.parent=parent
        self.dformat=self.parent.dformat
        self.buf=self.read(1024) #Read the whole block first
        self.namebuf=self.read(1024) #Read the name record next
        self.NEXT=int(struct.unpack(self.dformat,self.buf[ 0: 0+8])[0]) #Address of next summary record in file (zero if last)
        self.PREV=int(struct.unpack(self.dformat,self.buf[ 8: 8+8])[0]) #Address of prev summary record in file (zero if first)
        self.NSUM=int(struct.unpack(self.dformat,self.buf[16:16+8])[0]) #Number of summaries actually in this summary record
        self.SS=self.parent.SS
    def read(self,len):
        return self.parent.inf.read(len)
    def tell(self):
        return self.parent.inf.tell()
    def __iter__(self):
        for i in range(self.NSUM):
            result = daf_summary.make(self,self.parent,self.buf[8*(3+i*self.SS):8*(3+(i+1)*self.SS)],str(self.namebuf[i*self.parent.NC:(i+1)*self.parent.NC],encoding="cp437").strip(" \t\n\r\0"))
            yield result
    def __str__(self):
        result=(("NEXT: %d\n" +
                 "PREV: %s\n" +
                 "NSUM: %d") %
                (self.NEXT,self.PREV,self.NSUM))
        return result

class double_array_file:
    def __init__(self,inf):
        """
        Read a double-array file. This class reads the DAF header itself, and delegates reading the
        records, comments, and data. It therefore has fields representing the header, and lists
        of other class objects representing the data.
        :param inf: string filename or file-like object
        """
        try:
            self.inf=open(inf,"rb")
            self.needToClose=True
        except TypeError:
            self.inf=inf
            self.needToClose=False
        buf=self.inf.read(1024) #Read the whole header record first
        self.LOCIDW=str(buf[0:0+8],encoding='ASCII')
        if self.LOCIDW[0:4]!="DAF/":
            raise ValueError("Not a DAF: Header is %s"%self.header)
        self.subtype=self.LOCIDW[4:8].strip(" \t\n\r\0")
        self.LOCFMT=str(buf[88:88+8],encoding='ASCII') #IEEE format string
        if self.LOCFMT=="LTL-IEEE":
            self.endian="<"
        elif self.LOCFMT=="BIG-IEEE":
            self.endian=">"
        else:
            raise ValueError("Unknown double-precision format %s"%self.LOCFMT)
        self.iformat=self.endian+"i"
        self.dformat=self.endian+"d"
        self.ND=struct.unpack(self.iformat,buf[ 8: 8+4])[0] #Address  8, number of components in each array summary
        self.NI=struct.unpack(self.iformat,buf[12:12+4])[0] #Address 12, number of integer components in each array summary
        self.LOCIFN=str(buf[16:16+60],encoding='ASCII')     #Address 16, internal name or description of array file
        self.FWARD=struct.unpack(self.iformat,buf[76:76+4])[0] #Address 76, record number of initial summary in the file
        self.BWARD=struct.unpack(self.iformat,buf[80:80+4])[0] #Address 80, record number of final summary in the file
        self.FREE=struct.unpack(self.iformat,buf[84:84+4])[0]  #Address 84, first free address in file
        self.FTPSTR=str(buf[700:700+28],encoding='cp437') #Address 700, FTP validation string
        self.NC = 8 * (self.ND + (self.NI + 1) // 2)
        self.SS = self.ND + (self.NI + 1) // 2  # Size of a summary
        self.NS = 125 // self.SS  # Maximum number of summaries per summary record
    def __enter__(self):
        return self
    def __exit__(self,type,value,traceback):
        if self.needToClose:
            self.inf.close()
    def __str__(self):
        result=(("LOCIDW: %s\n"+
                 "ND:     %d\n"+
                 "NI:     %d\n" +
                 "NC*:    %d\n"+
                 "SS*:    %d\n" +
                 "NS*:    %d\n"+
                 "LOCIFN: %s\n" +
                 "FWARD:  %d\n" +
                 "BWARD:  %d\n" +
                 "FREE:   %d\n" +
                 "LOCFMT: %s") %
                (self.LOCIDW,self.ND,self.NI,self.NC,self.SS,self.NS,self.LOCIFN,self.FWARD,self.BWARD,self.FREE,self.LOCFMT))
        return result
    def __iter__(self):
        #Seek to the first summary record
        self.inf.seek((self.FWARD-1)*1024) #One for one-based, one for header record at beginning of file
        done=False
        while not done:
            result=daf_summary_record(self)
            yield result
            if result.NEXT==0:
                done=True
            else:
                self.inf.seek((result.NEXT-1)*1024)
    def dump(self,dump_lines=True):
        print(str(in_daf))
        for i,sr in enumerate(in_daf):
            print("Summary record %d"%i)
            print(str(sr))
            for j,sum in enumerate(sr):
                print("Summary %d" % j)
                print(str(sum))
                print(str(sum.segment()))
                print("i*,"+sum.segment().csvheader())
                if dump_lines:
                    for k,line in enumerate(sum.segment()):
                        print(k,","+str(line))

if __name__=="__main__":
    dump_lines=True
    from spiceypy import furnsh as cspice_furnsh
    cspice_furnsh("../../Data/spice/phoenix/phx.tsc")
    cspice_furnsh("../../Data/spice/phoenix/naif.tls")
    with double_array_file("../../Data/spice/phoenix/phx_edl_rec_att.bc") as in_daf:
        in_daf.dump(dump_lines)
    with double_array_file("../../Data/spice/insight/insight_nom_2016e09o_edl_v1.bsp") as in_daf:
        in_daf.dump(dump_lines)
    with double_array_file("../../Data/spice/phoenix/phx_edl_rec_traj.bsp") as in_daf:
        in_daf.dump(dump_lines)
