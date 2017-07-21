'''
Created on Jul 21, 2017

@author: chrisj
'''
import NetlistTree

class SN74xx_4x2gate(NetlistTree.Component):
    @staticmethod
    def factory(ref,comp):
        part=comp["part"]
        if   part=="74xx:74LS86"     : 
            return  SN7486         (ref,comp)
        elif part=="74xx:74LS08"     : 
            return  SN7408         (ref,comp)
        elif part=="74xx:74LS32"     : 
            return  SN7432         (ref,comp)
        else:
            print("No model found for "+part)
    def __init__(self,ref,comp):
        self.IN=[[False]*5,[False]*5]
        #          (0) 1  2  3  4  5  6  7  8  9 10 11 12 13 14
        self.unit=[-8, 1, 1,-1, 2, 2,-2,-9,-3, 3, 3,-4, 4, 4,-9]
        self.pin= [-1, 1, 2, 3, 1, 2, 3,-9, 3, 1, 2, 3, 1, 2,-9]
        self.out= [-1, 3, 6, 8,11]
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        #print("ref:",self.ref,"74xx:",self.value)
        #print("  IN-:",self.IN)
        pinnum=int(pin)
        unit=self.unit[pinnum]
        pintype=self.pin[pinnum]
        #print("  pinnum:",pinnum,"unit:",unit,"pintype:",pintype)
        if pintype==1 or pintype==2:
            #Input A or B
            if self.IN[pintype-1][unit]==value: return
            self.IN[pintype-1][unit]=value
            #print("  IN+:",self.IN)
        else:
            print("Trying to set an output or power")
        a=self.IN[0][unit]
        b=self.IN[1][unit]
        #print("  a:",a,"b:",b)
        newval=self.logic(a,b)
        newtime=time+self.propDelay
        outpin=str(self.out[unit])
        outnet=self.outputs[outpin]
        #print("  logic:",newval)
        #print("  eventtime:",newtime,"outpin:",outpin,"outnet:",outnet)
        context.addEvent(newtime,outnet,newval)

class SN7486(SN74xx_4x2gate):
    """
    Quad 2-input XOR
    """
    def logic(self,A,B):
        return A ^ B

class SN7408(SN74xx_4x2gate):
    """
    Quad 2-input AND
    """
    def logic(self,A,B):
        return A and B

class SN7432(SN74xx_4x2gate):
    """
    Quad 2-input OR
    """
    def logic(self,A,B):
        return A or B

