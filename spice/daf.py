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
        result ="%s\n"%self.name
        result+="SCLK0:    %20.15e" %self.sclk0
        try:
            et=cspice_sct2e(self.target//1000,self.sclk0)
            scdecd=cspice_scdecd(self.target//1000,self.sclk0)
            result+=" (%s ET=%16.6f %s)"%(scdecd,et,cspice_etcal(et))
        except:
            pass #If there is an error, it means that the correct kernels aren't available.
        result+="\n"
        result+="SCLK1:    %20.15e" %self.sclk1
        try:
            et=cspice_sct2e(self.target//1000,self.sclk1)
            scdecd=cspice_scdecd(self.target//1000,self.sclk1)
            result+=" (%s ET=%16.6f %s)"%(scdecd,et,cspice_etcal(et))
        except:
            pass #If there is an error, it means that the correct kernels aren't available.
        result+="\n"
        result+=(("Inst:  %7d\n"+
                  "Frame: %7d\n"+
                  "Type:  %7d\n"+
                  "Rates? %7d\n"+
                  "Addr0: %7d\n"+
                  "Addr1: %7d") %
                (self.target,self.frame,self.type,self.rates,self.addr0,self.addr1))
        return result
    def segment(self):
        return daf_CKSegment(self,self.daf)

class daf_SPKSegment:
    csvheaders=[0]*14
    csvheaders[1]="dt*,ET,ETCAL*,TL"
    for i in range(15):
        csvheaders[1]+=",G[%02d]"%(i)
    csvheaders[1]+=",X,DXDT,Y,DYDT,Z,DZDT"
    for i in range(3):
        for j in range(15):
            csvheaders[1]+=",DT[%02d%02d]"%(i,j)
    csvheaders[1]+=",KQMAX1,KQ[0],KQ[1],KQ[2]"
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
            return daf_SPKStateline(self.N,i,self.buf)
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
    csvheaders[3]="SCLK,SCLK,ET,ETCAL,q0,q1,q2,q3,wx,wy,wz"
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
            return daf_CK03line(self.summary.target//1000,self.N,i,self.buf,self.summary.rates)
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


class daf_SPKStateline:
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
        import numpy as np
        self.FD=np.zeros(15)*float('nan')

        # initialize arrays for spke01
        self.G = np.zeros(15)
        self.REFPOS = np.zeros(3)
        self.REFVEL = np.zeros(3)
        self.KQ = np.array([0, 0, 0])
        self.FC = np.zeros(15)
        self.FC[0] = 1.0
        self.WC = np.zeros(13)
        self.W = np.zeros(17)

        self.spke01(self.e[0],self.e)
#        self.eval(self.TL)
    def __str__(self):
        result="%30.14f,%s" %(self.et,cspice_etcal(self.et))
        for i in range(self.n_elements):
            result+=",%25.15e"%self.e[i]
        return result
    """
    def getTL(self): return self.e[0]
    def getG (self): return self.e[1:15+1]
    def getREFPOS(self): return self.e[16:22:2]
    def getREFVEL(self): return self.e[17:23:2]
    def getDT    (self): return self.e[22:22+45]
    def getKQMAX1(self): return int(self.e[67])
    def getKQ    (self): return self.e[68:71]
    TL    =property(fget=getTL    ,doc="Ephemeris time of end of segment")
    G     =property(fget=getG     ,doc="Stepsize function vector")
    REFPOS=property(fget=getREFPOS,doc="Reference position vector")
    REFVEL=property(fget=getREFVEL,doc="Reference velocity vector")
    DT    =property(fget=getDT    ,doc="Modified divided difference arrays")
    KQMAX1=property(fget=getKQMAX1,doc="Maximum integration order plus 1")
    KQ    =property(fget=getKQ    ,doc="Integration order array")
    """
    def getFC(self,J):
        print("Accessing element %d, value %f"%(J,self.FD[J]))
        return self.FD[J]
    def setFC(self,J,val):
        print("Writing element %d, was value %f, will be value %f"%(J,self.FD[J],val))
        self.FD[J]=val
    def eval(self,ET):
        import numpy as np
        STATE=np.zeros(6)
        DELTA = ET - self.TL
        TP = DELTA
        MQ2 = self.KQMAX1 - 2
        KS = self.KQMAX1 - 1
        self.setFC(0,1.0)
        WC = np.zeros(13)
        W  = np.zeros(17)

        #     This is clearly collecting some kind of coefficients.
        #     The problem is that we have no idea what they are...
        #
        #     The G coefficients are supposed to be some kind of step size
        #     vector.
        #
        #     TP starts out as the delta t between the request time
        #     and the time for which we last had a state in the MDL file.
        #     We then change it from DELTA  by the components of the stepsize
        #     vector G.
        for J in range(1, MQ2 + 1):
            self.setFC(J,TP / self.G[J-1])
            WC[J-1] = DELTA / self.G[J-1]
            TP = DELTA + self.G[J-1]
        #     Collect KQMAX1 reciprocals.
        for J in range(1, self.KQMAX1 + 1):
            W[J-1] = 1.0 / float(J)

        #     Compute the W(K) terms needed for the position interpolation
        #     (Note,  it is assumed throughout this routine that KS, which
        #     starts out as KQMAX1-1 (the ``maximum integration'')
        #     is at least 2.
        JX = 0
        KS1 = KS - 1
        while KS >= 2:
            JX = JX + 1
            for J in range(1, JX + 1):
                W[J + KS - 1] = self.getFC(J) * W[J + KS1 - 1] - WC[J - 1] * W[J + KS - 1]
            KS = KS1
            KS1 = KS1 - 1

        #     Perform position interpolation: (Note that KS = 1 right now.
        #     We don't know much more than that.)
        for I in range(3):
            KQQ = int(self.KQ[I])
            SUM = 0.0
            for J in range(KQQ, 0, -1):
                SUM = SUM + self.DT[J - 1+I*15] * W[J + KS - 1]
            STATE[I] = self.REFPOS[I] + DELTA * (self.REFVEL[I] + DELTA * SUM)

        #     Again we need to compute the W(K) coefficients that are
        #     going to be used in the velocity interpolation.
        #     (Note, at this point, KS = 1, KS1 = 0.)
        for J in range(1, JX + 1):
            W[J + KS - 1] = self.getFC(J) * W[J + KS1 - 1] - WC[J - 1] * W[J + KS - 1]
        KS = KS - 1

        #     Perform velocity interpolation:
        for I in range(1, 3 + 1):
            KQQ = int(self.KQ[I])
            SUM = 0.0

            for J in range(KQQ, 0, -1):
                SUM = SUM + self.DT[J - 1+I*15] * W[J + KS - 1]

            STATE[I + 3] = self.REFVEL[I] + DELTA * SUM

        #     That's all folks.  We don't know why we did anything, but
        #     at least we can tell structurally what we did.
        return STATE
    def test(self,et0,N=100):
        """
        Compare the calculation of the state in pure python to that provided by spiceypy (which is an interface
        to the cspice library). This should be near bit-perfect.
        :param ET0: Time to begin at. Since the record only contains the end time, we need a beginning time. The
                    beginning time *might* be encoded in the DeltaT's, but this isn't known for sure yet
        :param N: Number of points to evaluate at. If you pass 100, you will get the point at the beginning,
                  99 points in between, and the point at the end, so a total of 101 points.
        :return: Sum of squares of differences of position components. Ideally should be zero, should be only a couple
                 of ulps.
        """
        from spiceypy import spkezr as cspice_spkezr
        import numpy as np
        et1=self.et
        result=0.0
        for i in range(N):
            et=et0+i*(et1-et0)/N
            #et=et1
            ref=cspice_spkezr("-189",et,"J2000","NONE","399")[0][0:3]
            test=np.array(self.spke01(et,self.e)[0:3])
            result+=np.sum((ref-test)**2)
        return result
    def spke01(self, ET, RECORD):
        import numpy as np
        """Compute position and velocity from a Modified Difference Array record

        Inputs:
            ET: Epoch time to evaluate position and velocity (seconds since J2000)
            RECORD: A record of Modified Difference Array
        Returns: STATE
            STATE: A numpy array which contains position and velocity
        """

        # This method has been translated from SPKE01 of SPICE Toolkit and
        # modified by Shushi Uetsuki.
        #
        # SPICE Toolkit for FORTRAN : http://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
        # SPK Required Reading : http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
        #
        # Original FORTRAN code uses 'SAVE' directive, and it means all variable
        # should be saved for next call.  So i decided to make almost all
        # variables to be instance variable.  Some of them are initialized in
        # __init__ method.

        STATE = np.zeros(6)

        #     Variable  I/O  Description
        #     --------  ---  --------------------------------------------------
        #     ET         I   Target epoch.
        #     RECORD     I   Data record.
        #     STATE      O   State (position and velocity).
        #
        # $ Detailed_Input
        #
        #     ET          is a target epoch, at which a state vector is to
        #                 be computed.
        #
        #     RECORD      is a data record which, when evaluated at epoch ET,
        #                 will give the state (position and velocity) of some
        #                 body, relative to some center, in some inertial
        #                 reference frame.
        #
        # $ Detailed_Output
        #
        #     STATE       is the state. Units are km and km/sec.
        #
        # $ Parameters
        #
        #     None.
        #
        # $ Exceptions
        #
        #     None.
        #
        # $ Files
        #
        #     None.
        #
        # $ Particulars
        #
        #     The exact format and structure of type 1 (difference lines)
        #     segments are described in the SPK Required Reading file.
        #
        #     Difference lines (DL's) are generated by JPL navigation
        #     system programs P and PV. Each data record is equivalent
        #     to the (slightly rearranged) 'P' portion of a NAVIO PV file
        #     data record.
        #
        #     SPKE01 is a specialized version of Fred Krogh's subroutine DAINT.
        #     Only the calling sequence has been changed.
        #
        #     Because the original version was undocumented, only Fred
        #     knows how this really works.
        #
        # $ Examples
        #
        #     None.
        #
        # $ Restrictions
        #
        #     Unknown.
        #
        # $ Literature_References
        #
        #     NAIF Document 168.0, "S- and P- Kernel (SPK) Specification and
        #     User's Guide"
        #
        # $ Author_and_Institution
        #
        #     F.T. Krogh      (JPL)
        #     I.M. Underwood  (JPL)
        #
        # $ Version
        #
        # -    SPICELIB Version 1.1.0, 14-FEB-1997 (WLT)
        #
        #        The goto's were removed and loop and if structures
        #        revealed.  We still don't know exactly what's going
        #        on, but at least the bones of this routine have been
        #        cleaned off and are ready for assembly. (WLT)
        #
        # -    SPICELIB Version 1.0.4, 30-OCT-1996 (WLT)
        #
        #        Removed redundant SAVE statements from the declaration
        #        section.  Thanks to Steve Schlaifer for finding this
        #        error.
        #
        # -    SPICELIB Version 1.0.3, 10-MAR-1992 (WLT)
        #
        #        Comment section for permuted index source lines was added
        #        following the header.
        #
        # -    SPICELIB Version 1.0.2, 23-AUG-1991 (HAN)
        #
        #        SPK01 was removed from the Required_Reading section of the
        #        header. The information in the SPK01 Required Reading file
        #        is now part of the SPK Required Reading file.
        #
        # -    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN)
        #
        #        Literature references added to the header.
        #
        # -    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) (FTK)
        #
        # -&
        #
        # $ Index_Entries
        #
        #     evaluate type_1 spk segment
        #
        # -&



        #
        #     Unpack the contents of the MDA array.
        #
        #        Name    Dimension  Description
        #        ------  ---------  -------------------------------
        #        TL              1  Final epoch of record
        #        G              15  Stepsize function vector
        #        REFPOS          3  Reference position vector
        #        REFVEL          3  Reference velocity vector
        #        DT         15,NTE  Modified divided difference arrays
        #        KQMAX1          1  Maximum integration order plus 1
        #        KQ            NTE  Integration order array
        #
        #     For our purposes, NTE is always 3.
        #
        self.TL = RECORD[0]
        self.G = RECORD[1:1 + 15]
        #
        #     Collect the reference position and velocity.
        #
        self.REFPOS[0] = RECORD[16]
        self.REFVEL[0] = RECORD[17]

        self.REFPOS[1] = RECORD[18]
        self.REFVEL[1] = RECORD[19]

        self.REFPOS[2] = RECORD[20]
        self.REFVEL[2] = RECORD[21]

        self.DT = np.reshape(RECORD[22:22 + 45], (15, 3), order='F')

        self.KQMAX1 = int(RECORD[67])
        self.KQ[0] = int(RECORD[68])
        self.KQ[1] = int(RECORD[69])
        self.KQ[2] = int(RECORD[70])
        #
        #     Next we set up for the computation of the various differences
        #
        self.DELTA = ET - self.TL
        self.TP = self.DELTA
        self.MQ2 = self.KQMAX1 - 2
        self.KS = self.KQMAX1 - 1
        #
        #     This is clearly collecting some kind of coefficients.
        #     The problem is that we have no idea what they are...
        #
        #     The G coefficients are supposed to be some kind of step size
        #     vector.
        #
        #     TP starts out as the delta t between the request time
        #     and the time for which we last had a state in the MDL file.
        #     We then change it from DELTA  by the components of the stepsize
        #     vector G.
        #
        for J in range(1, self.MQ2 + 1):
            self.FC[J] = self.TP / self.G[J - 1]
            self.WC[J - 1] = self.DELTA / self.G[J - 1]
            self.TP = self.DELTA + self.G[J - 1]
        #
        #     Collect KQMAX1 reciprocals.
        #
        for J in range(1, self.KQMAX1 + 1):
            self.W[J - 1] = 1.0 / float(J)
        #
        #     Compute the W(K) terms needed for the position interpolation
        #     (Note,  it is assumed throughout this routine that KS, which
        #     starts out as KQMAX1-1 (the ``maximum integration'')
        #     is at least 2.
        #
        self.JX = 0
        self.KS1 = self.KS - 1

        while self.KS >= 2:

            self.JX = self.JX + 1

            for J in range(1, self.JX + 1):
                self.W[J + self.KS - 1] = self.FC[J] * self.W[J + self.KS1 - 1] - self.WC[J - 1] * self.W[
                    J + self.KS - 1]

            self.KS = self.KS1
            self.KS1 = self.KS1 - 1
        #
        #     Perform position interpolation: (Note that KS = 1 right now.
        #     We don't know much more than that.)
        #
        for I in range(1, 3 + 1):

            self.KQQ = self.KQ[I - 1]
            self.SUM = 0.0

            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J - 1, I - 1] * self.W[J + self.KS - 1]

            #x=x0+delta*dx0+delta**2*sum
            STATE[I - 1] = self.REFPOS[I - 1] + self.DELTA * (self.REFVEL[I - 1] + self.DELTA * self.SUM)
        #
        #     Again we need to compute the W(K) coefficients that are
        #     going to be used in the velocity interpolation.
        #     (Note, at this point, KS = 1, KS1 = 0.)
        #
        for J in range(1, self.JX + 1):
            self.W[J + self.KS - 1] = self.FC[J] * self.W[J + self.KS1 - 1] - self.WC[J - 1] * self.W[J + self.KS - 1]

        self.KS = self.KS - 1

        #
        #     Perform velocity interpolation:
        #
        for I in range(1, 3 + 1):
            self.KQQ = self.KQ[I - 1]
            self.SUM = 0.0

            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J - 1, I - 1] * self.W[J + self.KS - 1]

            STATE[I + 3 - 1] = self.REFVEL[I - 1] + self.DELTA * self.SUM

        #
        #     That's all folks.  We don't know why we did anything, but
        #     at least we can tell structurally what we did.
        #

        return STATE

class daf_CK03line:
    def __init__(self,clknum,N,i,buf,rates):
        self.clknum=clknum
        self.rates=rates
        if self.rates:
            self.n_elements = 7
            self.w = buf[i * self.n_elements+4:(i + 1) * self.n_elements]
        else:
            self.n_elements=4
        self.q = buf[i * self.n_elements:i * self.n_elements+4]

        self.sclk=buf[N*self.n_elements+i]
    def __str__(self):
        result="%25.15e"%self.sclk
        try:
            et=cspice_sct2e(self.clknum, self.sclk)
            result+=",%s,%16.6f,%s" %(cspice_scdecd(self.clknum,self.sclk),et,cspice_etcal(et))
        except:
            result+=",,,"
        for i in range(4):
            result+=",%25.15e"%self.q[i]
        if self.rates:
            for i in range(3):
                result += ",%25.15e" % self.w[i]
        else:
            result+=",,,"
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
    def dump(self,dump_lines=True,dump_comments=True):
        print(str(self))
        if dump_comments:
            for line in self.comments():
                print(line)
        for i,sr in enumerate(self):
            print("Summary record %d"%i)
            print(str(sr))
            for j,sum in enumerate(sr):
                print("Summary %d" % j)
                print(str(sum))
                print(str(sum.segment()))
                print("i*,"+sum.segment().csvheader())
                if dump_lines:
                    et0=sum.et0
                    for k,line in enumerate(sum.segment()):
                        dt=line.et-et0
                        print(k,",%15.6f,"%dt,str(line))
#                        tval=line.test(et0)
#                        if True: #(tval>1e-10):
#                            print(line.et-et0,",",tval)
                        et0=line.et
    def comments(self):
        block=1
        done=False
        line=""
        while not done:
            self.inf.seek(block*1024)
            bb=self.inf.read(1000)
            for b in bb:
                if b==0:
                    yield line
                    line=""
                elif b==4:
                    yield line
                    return
                else:
                    line+=chr(b)
            block+=1

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description="Translate a Spice Double-Array File (DAF) into plain text")
    parser.add_argument('-x','--records',action='store_true',help="Dump all data records. May result in a very large amount of output!")
    parser.add_argument('-r','--comment',action='store_true',help="Dump file comments.")
    parser.add_argument('-a','--all',action='store_true',help="Dump everything that can be dumped")
    parser.add_argument('inf',nargs=1,help="Kernel to translate. Must be a bc or bsp kernel.")
    parser.add_argument('kernel',nargs='*',help="Kernels to furnish. Typically these are sclk and/or lsk files, needed to interpret bc timestamps")
    result=parser.parse_args()
    dump_lines=result.records or result.all
    dump_comments=result.comment or result.all
    inf=result.inf[0]

    from spiceypy import furnsh as cspice_furnsh
    for furnsh in result.kernel:
        cspice_furnsh(furnsh)
    with double_array_file(inf) as in_daf:
        cspice_furnsh(inf) #Just to run the test
        in_daf.dump(dump_lines=dump_lines,dump_comments=dump_comments)
#    with double_array_file("../../Data/spice/insight/insight_nom_2016e09o_edl_v1.bsp") as in_daf:
#        in_daf.dump(dump_lines)
#    with double_array_file("../../Data/spice/phoenix/phx_edl_rec_traj.bsp") as in_daf:
#        in_daf.dump(dump_lines)
