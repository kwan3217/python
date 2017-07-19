import antlr4
import heapq
from Netlist.NetlistLexer   import NetlistLexer 
from Netlist.NetlistParser  import NetlistParser 
from Netlist.NetlistVisitor import NetlistVisitor 
 
class NetlistTree(NetlistVisitor):
    @staticmethod
    def sep(children):
        dicts={}
        nondicts=[]
        for part in children:
            if isinstance(part,dict):
                dicts.update(part)
            else:
                nondicts.append(part)
        return (dicts,nondicts)        
    def visitTerminal(self, node):
        return node.symbol.text
    def aggregateResult(self, aggregate, nextResult):
        if isinstance(aggregate,dict):
            if isinstance(nextResult,dict):
                aggregate.update(nextResult)
        elif isinstance(aggregate,list):
            aggregate.append(nextResult)
        else:
            aggregate=[nextResult]
        return aggregate
    def visitString(self, ctx:NetlistParser.StringContext):
        result=ctx.children[0].symbol.text
        if result[0]=='"':
            result=result[1:-1]
        return result
    def visitBoringNode(self, ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:result[1]}
    def visitBoringIntNode(self, ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:int(result[1])}
    def visitOptionalBoringNode(self, ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        if len(result) < 3:
            return {result[0][1:]:""}
        return {result[0][1:]:result[1]}
    #Nontrivial Nonterminals
    def visitNetlist(self,ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        return result[0]
    def visitExport(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {'components':dicts["components"],
                'libparts'  :dicts["libparts"  ],
                'nets'      :dicts["nets"      ]}
    def visitComponents(self,ctx:NetlistParser.StringContext):
        children=self.visitChildren(ctx)
        components={}
        for component in children[1:-1]:
            ref=component["ref"]
            del component["ref"]
            components[ref]=component
        return {"components":components}
    def visitComp(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitLibsource(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libsource":dicts,"part":dicts["lib"]+":"+dicts["part"]}
    def visitLibparts(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libparts":dicts}
    def visitLibpart(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {dicts["lib"]+":"+dicts["part"]:dicts}
    def visitPins(self,ctx:NetlistParser.StringContext):
        children=self.visitChildren(ctx)
        pins={}
        for pin in children[1:-1]:
            num=pin["num"]
            del pin["num"]
            pins[num]=pin
        return {"pins":pins}
    def visitPin(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitFields(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"fields":dicts}
    def visitField(self,ctx:NetlistParser.StringContext):
        (dicts,nondicts)=self.sep(self.visitChildren(ctx))
        return {dicts["name"]:nondicts[1]}
    def visitNets(self,ctx:NetlistParser.StringContext):
        children=self.visitChildren(ctx)
        nets={}
        for net in children[1:-1]:
            name=net["name"]
            del net["name"]
            nets[name]=net
        return {"nets":nets}
    def visitNet(self,ctx:NetlistParser.StringContext):
        (dicts,nondicts)=self.sep(self.visitChildren(ctx))
        nodes=[]
        for node in nondicts[1:-1]:
            (nodedict,_)=self.sep(node)
            nodes.append(nodedict)
        dicts.update({"nodes":nodes})
        return dicts
    
    #Don't care nonterminals
    def visitDesign(self,ctx:NetlistParser.StringContext): return None
    #Boring nodes
    def visitAlias      (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitCode       (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitDescription(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitDocs       (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitFootprint  (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitFp         (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitLib        (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitLogical    (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitName       (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitNames      (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitNum        (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitNumber     (self,ctx:NetlistParser.StringContext):return self.visitBoringIntNode(ctx)
    def visitPart       (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitPin__innode(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitRef        (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitSource     (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitTool       (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamp     (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamps    (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitType_      (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitUri        (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitValue      (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitVersion    (self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    #Optional boring nodes
    def visitCompany(self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitDate   (self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitRev    (self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitTitle  (self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)

class Component:
    @staticmethod
    def factory(ref,comp):
        part=comp["part"]
        if part=="kwire:AND"      : return KwireAND      (ref,comp)
        if part=="kwire:OR"       : return KwireOR       (ref,comp)
        if part=="kwire:NOT"      : return KwireNOT      (ref,comp)
        if part=="kwire:XOR"      : return KwireXOR      (ref,comp)
        if part=="kwire:CONTROL"  : return KwireCONTROL  (ref,comp)
        if part=="kwire:INDICATOR": return KwireINDICATOR(ref,comp)
    def __init__(self,ref,comp):
        self.ref=ref
        self.libsource=comp["libsource"]
        self.part=comp["part"]
        self.value=comp["value"]
        self.outputs={}
        self.propDelay=0
    def addOutput(self,name,net):
        self.outputs[name]=net
    def setInput(self,context,time,pin,value):
        pass
    def __str__(self):
        return self.__dict__.__str__()
    
class KwireAND(Component):
    def __init__(self,ref,comp):
        self.A=False
        self.B=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.A: return
            self.A=value
        else:
            if value==self.B: return
            self.B=value
        context.addEvent(time+self.propDelay,self.outputs["3"],self.A and self.B)

class KwireOR(Component):
    def __init__(self,ref,comp):
        self.A=False
        self.B=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.A: return
            self.A=value
        else:
            if value==self.B: return
            self.B=value
        context.addEvent(time+self.propDelay,self.outputs["3"],self.A or self.B)

class KwireNOT(Component):
    def __init__(self,ref,comp):
        self.A=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.A: return
        self.A=value
        context.addEvent(time+self.propDelay,self.outputs["2"],not self.A)

class KwireXOR(Component):
    def __init__(self,ref,comp):
        self.A=False
        self.B=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if pin=='1':
            if value==self.A: return
            self.A=value
        else:
            if value==self.B: return
            self.B=value
        context.addEvent(time+self.propDelay,self.outputs["3"],self.A ^ self.B)

class KwireCONTROL(Component):
    def __init__(self,ref,comp):
        self.Y=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.Y: return
        self.Y=value
        context.addEvent(time+self.propDelay,self.outputs["1"],value)

class KwireINDICATOR(Component):
    def __init__(self,ref,comp):
        self.A=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.A: return
        self.A=value
        
class Event:
    def __init__(self,time,net,value):
        self.time=time
        self.net=net
        self.value=value
    def __lt__(self,other):
        return self.time<other.time

class Context:
    def __init__(self,infn):
        inf = antlr4.FileStream(infn)
        lexer = NetlistLexer(inf)
        stream = antlr4.CommonTokenStream(lexer)
        parser = NetlistParser(stream)
        tree = parser.netlist()
        visitor=NetlistTree()
        netlist=visitor.visit(tree)
        self.time=0
        self.eventList=[]
        self.nets=netlist["nets"]
        self.libparts=netlist["libparts"]
        self.components={}
        for ref,comp in netlist["components"].items():
            self.components[ref]=Component.factory(ref,comp)
        #Figure out the input for each net
        for name,net in self.nets.items():
            #For each node, figure out if the pin is an output
            for node in net["nodes"][:]: #The [:] makes a quick copy of the list in order to iterate on it
                comp=self.components[node["ref"]]
                libpart=self.libparts[comp.part]
                pin=libpart["pins"][node["pin"]]
                if pin["type"]=="output":
                    comp.addOutput(node["pin"],name)
                    net["nodes"].remove(node)
                    
    def dispatchEvent(self):
        if len(self.eventList)==0:
            return False
        #Find the event with the earliest time
        event=heapq.heappop(self.eventList)
        self.time=event.time
        net=event.net
        for node in self.nets[net]["nodes"]: 
            comp=self.components[node["ref"]]
            comp.setInput(context,context.time,node["pin"],event.value)
        return True
    def addEvent(self,time,net=None,value=False):
        #HiZ by definition does't affect the net, so they shouldn't 
        #create an event
        if value=='Z': return;
        #If there is another event with the
        #exact same name and time, then that 
        #event is replaced with this one.
        for event in self.eventList:
            if event.time==time and event.net==net:
                event.value=value
                return
        heapq.heappush(self.eventList,Event(time,net,value))
    
if __name__ == '__main__':
    context=Context('../../kicad/kwire/1bitdriver.net')
    inputs=[{"time": 1,"ref":"S101","value":True },
            {"time": 2,"ref":"S102","value":True },
            {"time": 3,"ref":"S101","value":False},
            {"time": 4,"ref":"S103","value":True },
            {"time": 5,"ref":"S101","value":True },
            {"time": 6,"ref":"S102","value":False},
            {"time": 7,"ref":"S101","value":False},
            {"time": 8,"ref":"S104","value":True },
            {"time": 9,"ref":"S101","value":True },
            {"time":10,"ref":"S102","value":True },
            {"time":11,"ref":"S101","value":False},
            {"time":12,"ref":"S103","value":False},
            {"time":13,"ref":"S101","value":True },
            {"time":14,"ref":"S102","value":False},
            {"time":15,"ref":"S101","value":False},
            {"time":16,"ref":"S104","value":False}]
    print(" T|- C B A|S K")
    print("--+-------+---")
    for inp in inputs:
        context.components[inp["ref"]].setInput(context,inp["time"],"1",inp["value"])
        while context.dispatchEvent():
#            for event in context.eventList:
#                print(event)
#                pass
            pass
        print("%2d|%d %d %d %d|%d %d" % (context.time,
              1 if context.components["S104"].Y else 0,
              1 if context.components["S103"].Y else 0,
              1 if context.components["S102"].Y else 0,
              1 if context.components["S101"].Y else 0,
              1 if context.components["D101"].A else 0,
              1 if context.components["D102"].A else 0))
