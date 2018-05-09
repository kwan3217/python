"""
Code for parsing Spice Double-Array Files (DAF). This is a superset of CK and SPK files. Note that while this
parses Spice files, it only makes minimal use of the spice libraries (date conversion only)
"""

import struct

class daf_summary:
    @staticmethod
    def make(sr,daf,buf,name):
        if daf.subtype=="SPK":
            return daf_SPKsummary(sr,daf,buf,name)
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
                 "ET0:    %30.14f\n"+
                 "ET1:    %30.14f\n"+
                 "Target: %d\n"+
                 "Center: %d\n"+
                 "Frame:  %d\n"+
                 "Type:   %d\n"+
                 "Addr0:  %d\n"+
                 "Addr1:  %d") %
                (self.name,self.et0,self.et1,self.target,self.center,self.frame,self.type,self.addr0,self.addr1))
        return result
    def segment(self):
        if self.type==9:
            return daf_SPKSegment(self,self.daf)

class daf_SPKSegment:
    csvheaders=[0]*10
    csvheaders[9]="ET,x,y,z,dxdt,dydt,dzdt"
    def __init__(self,summary,daf):
        self.summary=summary
        self.daf=daf
        daf.inf.seek(8*(summary.addr0-1))
        segment_size=summary.addr1-summary.addr0+1
        format="%s%d%s"%(daf.dformat[0],segment_size,daf.dformat[1])
        self.buf=struct.unpack(format,daf.inf.read(segment_size*8))
        self.N=int(self.buf[-1])
    def line(self,i):
        if self.summary.type==9:
            return daf_SPK09line(self.N,i,self.buf)
    def __iter__(self):
        for i in range(self.N):
            result = self.line(i)
            yield result
    def csvheader(self):
        return self.csvheaders[self.summary.type]

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
        result="%30.14f,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e"%(self.et,self.x,self.y,self.z,self.dxdt,self.dydt,self.dzdt)
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
        self.dformat="<d"
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

if __name__=="__main__":
    with double_array_file("/home/chrisj/workspace/Data/spice/insight/insight_nom_2016e09o_edl_v1.bsp") as in_daf:
        print(str(in_daf))
        for i,sr in enumerate(in_daf):
            print("Summary record %d"%i)
            print(str(sr))
            for j,sum in enumerate(sr):
                print("Summary %d" % j)
                print(str(sum))
                print("i,"+sum.segment().csvheader())
                for k,line in enumerate(sum.segment()):
                    print(k,","+str(line))
