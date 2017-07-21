'''
Created on Jul 21, 2017

@author: chrisj
'''

import NetlistTree

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
        elif pin=='2':
            if value==self.S: return
            self.S=value
        else:
            print("Trying to write to transistor drain")
        if self.G==True:
            result=self.S
        else:
            result='Z'
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
        elif pin=='2':
            if value==self.S: return
            self.S=value
        else:
            print("Trying to write to transistor drain")
        if self.G==False:
            result=self.S
        else:
            result='Z'
        context.addEvent(time+self.propDelay,self.outputs["3"],result)

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
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.Y[17-int(pin)]: return
        self.Y[17-int(pin)]=value
        context.addEvent(time+self.propDelay,self.outputs[pin],value)
        #print(self.ref,self.Y)

class Kwire8INDICATOR(KwireComponent):
    def __init__(self,ref,comp):
        self.A=['Z']*(8+1)
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.A[17-int(pin)]: return
        self.A[17-int(pin)]=value
        #print(self.ref,self.A)


