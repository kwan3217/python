'''
Created on Jul 21, 2017

@author: chrisj
'''

import NetlistTree
from gi.overrides.keysyms import value

class KwireComponent(NetlistTree.Component):
    @staticmethod
    def factory(ref,comp):
        part=comp["part"]
        if   part=="kwire:AND"       :
            return  KwireAND       (ref,comp)
        elif part=="kwire:OR"        : 
            return  KwireOR        (ref,comp)
        elif part=="kwire:NOT"       : 
            return  KwireNOT       (ref,comp)
        elif part=="kwire:XOR"       : 
            return  KwireXOR       (ref,comp)
        elif part=="kwire:CONTROL"   : 
            return  KwireCONTROL   (ref,comp)
        elif part=="kwire:8CONTROL"  : 
            return  Kwire8CONTROL  (ref,comp)
        elif part=="kwire:8INDICATOR": 
            return  Kwire8INDICATOR(ref,comp)
        elif part=="kwire:NMOS":
            return  KwireNMOS(ref,comp)
        elif part=="kwire:PMOS":
            return  KwirePMOS(ref,comp)
        elif part=="kwire:NMOS_BODY":
            return  KwireNMOS_BODY(ref,comp)
        elif part=="kwire:PMOS_BODY":
            return  KwirePMOS_BODY(ref,comp)
        else:
            print("No model found for "+part)
            return KwireComponent(ref,comp)

class KwireAND(KwireComponent):
    def __init__(self,ref,comp):
        self.A='Z'
        self.B='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.A: return
            self.A=value
        else:
            if value==self.B: return
            self.B=value
        context.addEvent(time+self.propDelay,self.outputs["3"],self.A and self.B)

class KwireNMOS(KwireComponent):
    def __init__(self,ref,comp):
        self.S='Z'
        self.G='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.G: return
            self.G=value
            print(self.ref+" (NMOS) Setting gate to ",value)
        elif pin=='2':
            if value==self.S: return
            self.S=value
            print(self.ref+" (NMOS) Setting source to ",value)
        else:
            print("Trying to write to transistor drain")
        if self.G==True and self.S==False:
            result=self.S
            print(self.ref+" (NMOS) is on  (gate=",self.G,", source=",self.S,"), setting drain to ",result)
        else:
            result='Z'
            print(self.ref+" (NMOS) is off (gate=",self.G,", source=",self.S,"), setting drain to ",result)
        context.addEvent(time+self.propDelay,self.outputs["3"],result)

class KwireNMOS_BODY(KwireComponent):
    def __init__(self,ref,comp):
        self.S='Z'
        self.G='Z'
        self.B='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.B: return
            self.B=value
            print(self.ref+" (NMOS_BODY) Setting body to ",value)
        elif pin=='4':
            if value==self.G: return
            self.G=value
            print(self.ref+" (NMOS_BODY) Setting gate to ",value)
        elif pin=='2':
            if value==self.S: return
            self.S=value
            print(self.ref+" (NMOS_BODY) Setting source to ",value)
        else:
            print("Trying to write to transistor drain")
        if self.G==True and self.B==False:
            result=self.S
            print(self.ref+" (NMOS_BODY) is on  (gate=",self.G,", body=",self.B,"), setting drain to ",result)
        else:
            result='Z'
            print(self.ref+" (NMOS_BODY) is off (gate=",self.G,", body=",self.B,"), setting drain to ",result)
        context.addEvent(time+self.propDelay,self.outputs["3"],result)

class KwirePMOS(KwireComponent):
    def __init__(self,ref,comp):
        self.S='Z'
        self.G='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.G: return
            self.G=value
            print(self.ref+" (PMOS) Setting gate to ",value)
        elif pin=='2':
            if value==self.S: return
            self.S=value
            print(self.ref+" (PMOS) Setting source to ",value)
        else:
            print("Trying to write to transistor drain")
        if self.G==False and self.S==True:
            result=self.S
            print(self.ref+" (PMOS) is on  (gate=",self.G,", source=",self.S,"), setting drain to ",result)
        else:
            result='Z'
            print(self.ref+" (PMOS) is off (gate=",self.G,", source=",self.S,"), setting drain to ",result)
        context.addEvent(time+self.propDelay,self.outputs["3"],result)

class KwirePMOS_BODY(KwireComponent):
    def __init__(self,ref,comp):
        self.S='Z'
        self.G='Z'
        self.B='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.B: return
            self.B=value
            print(self.ref+" (PMOS) Setting body to ",value)
        elif pin=='3':
            if value==self.G: return
            self.G=value
            print(self.ref+" (PMOS) Setting gate to ",value)
        elif pin=='4':
            if value==self.S: return
            self.S=value
            print(self.ref+" (PMOS) Setting source to ",value)
        else:
            print("Trying to write to transistor drain")
        if self.G==False and self.B==True:
            result=self.S
            print(self.ref+" (PMOS_BODY) is on  (gate=",self.G,", body=",self.B,"), setting drain to ",result)
        else:
            result='Z'
            print(self.ref+" (PMOS_BODY) is off (gate=",self.G,", body=",self.B,"), setting drain to ",result)
        context.addEvent(time+self.propDelay,self.outputs["2"],result)
        
class KwireOR(KwireComponent):
    def __init__(self,ref,comp):
        self.A='Z'
        self.B='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.A: return
            self.A=value
        else:
            if value==self.B: return
            self.B=value
        context.addEvent(time+self.propDelay,self.outputs["3"],self.A or self.B)

class KwireNOT(KwireComponent):
    def __init__(self,ref,comp):
        self.A='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.A: return
        self.A=value
        context.addEvent(time+self.propDelay,self.outputs["2"],not self.A)

class KwireXOR(KwireComponent):
    def __init__(self,ref,comp):
        self.A='Z'
        self.B='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.A: return
            self.A=value
        else:
            if value==self.B: return
            self.B=value
        context.addEvent(time+self.propDelay,self.outputs["3"],self.A ^ self.B)

class KwireCONTROL(KwireComponent):
    def __init__(self,ref,comp):
        self.Y='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.Y: return
        self.Y=value
        context.addEvent(time+self.propDelay,self.outputs["1"],value)
        #print(self.ref,self.Y)

class Kwire8CONTROL(KwireComponent):
    def __init__(self,ref,comp):
        self.Y=['Z']*(8+1)
        self.pin1='Z'
        self.pin8='Z'
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        pinnum=int(pin)
        if pinnum<=8:
            if pinnum==1:
                print(self.ref,"Setting pin1 to ",value)
                self.pin1=value
            elif pinnum==8:
                print(self.ref,"Setting pin8 to ",value)
                self.pin8=value
            else:
                print("Input to a NC pin")
            return
            for i in range(8):
                swnum=i+1                
                if value=='Z':
                    outval='Z'
                elif value==True:
                    outval=self.pin1
                elif value==False:
                    outval=self.pin8
                context.addEvent(time+self.propDelay,self.outputs[pin],value)
        swnum=17-pinnum
        if value==self.Y[swnum]: return
        self.Y[swnum]=value
        if value=='Z':
            outval='Z'
        elif value==True:
            outval=self.pin1
        elif value==False:
            outval=self.pin8
        context.addEvent(time+self.propDelay,self.outputs[pin],value)
        #print(self.ref,self.Y)

class Kwire8INDICATOR(KwireComponent):
    def __init__(self,ref,comp):
        self.A=['Z']*(8+1)
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.A[17-int(pin)]: return
        self.A[17-int(pin)]=value
        print(self.ref,self.A)


