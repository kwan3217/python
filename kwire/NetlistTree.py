import antlr4
import heapq
from NetlistFile.NetlistFileLexer   import NetlistFileLexer 
from NetlistFile.NetlistFileParser  import NetlistFileParser
from NetlistFile.NetlistFileVisitor import NetlistFileVisitor
 
class NetlistTree(NetlistFileVisitor):
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
    def visitString(self, ctx:NetlistFileParser.StringContext):
        result=ctx.children[0].symbol.text
        if result[0]=='"':
            result=result[1:-1]
        return result
    def visitBoringNode(self, ctx:NetlistFileParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:result[1]}
    def visitBoringIntNode(self, ctx:NetlistFileParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:int(result[1])}
    def visitOptionalBoringNode(self, ctx:NetlistFileParser.StringContext):
        result=self.visitChildren(ctx)
        if len(result) < 3:
            return {result[0][1:]:""}
        return {result[0][1:]:result[1]}
    #Nontrivial Nonterminals
    def visitNetlist(self,ctx:NetlistFileParser.StringContext):
        result=self.visitChildren(ctx)
        return result[0]
    def visitExport(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {'components':dicts["components"],
                'libparts'  :dicts["libparts"  ],
                'nets'      :dicts["nets"      ]}
    def visitComponents(self,ctx:NetlistFileParser.StringContext):
        children=self.visitChildren(ctx)
        components={}
        for component in children[1:-1]:
            ref=component["ref"]
            del component["ref"]
            components[ref]=component
        return {"components":components}
    def visitComp(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitLibsource(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libsource":dicts,"part":dicts["lib"]+":"+dicts["part"]}
    def visitLibparts(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libparts":dicts}
    def visitLibpart(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {dicts["lib"]+":"+dicts["part"]:dicts}
    def visitPins(self,ctx:NetlistFileParser.StringContext):
        children=self.visitChildren(ctx)
        pins={}
        for pin in children[1:-1]:
            num=pin["num"]
            del pin["num"]
            pins[num]=pin
        return {"pins":pins}
    def visitPin(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitFields(self,ctx:NetlistFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"fields":dicts}
    def visitField(self,ctx:NetlistFileParser.StringContext):
        (dicts,nondicts)=self.sep(self.visitChildren(ctx))
        return {dicts["name"]:nondicts[1]}
    def visitNets(self,ctx:NetlistFileParser.StringContext):
        children=self.visitChildren(ctx)
        nets={}
        for net in children[1:-1]:
            name=net["name"]
            del net["name"]
            nets[name]=net
        return {"nets":nets}
    def visitNet(self,ctx:NetlistFileParser.StringContext):
        (dicts,nondicts)=self.sep(self.visitChildren(ctx))
        nodes=[]
        for node in nondicts[1:-1]:
            (nodedict,_)=self.sep(node)
            nodes.append(nodedict)
        dicts.update({"nodes":nodes})
        return dicts
    
    #Don't care nonterminals
    def visitDesign(self,ctx:NetlistFileParser.StringContext): return None
    #Boring nodes
    def visitAlias      (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitCode       (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitDescription(self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitDocs       (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitFootprint  (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitFp         (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitLib        (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitLogical    (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitName       (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitNames      (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitNum        (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitNumber     (self,ctx:NetlistFileParser.StringContext):return self.visitBoringIntNode(ctx)
    def visitPart       (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitPin__innode(self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitRef        (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitSource     (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitTool       (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamp     (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamps    (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitType_      (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitUri        (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitValue      (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitVersion    (self,ctx:NetlistFileParser.StringContext):return self.visitBoringNode(ctx)
    #Optional boring nodes
    def visitCompany(self,ctx:NetlistFileParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitDate   (self,ctx:NetlistFileParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitRev    (self,ctx:NetlistFileParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitTitle  (self,ctx:NetlistFileParser.StringContext):return self.visitOptionalBoringNode(ctx)

class Component:
    @staticmethod
    def factory(ref,comp):
        import SN74xx
        import KwireComponents
        part=comp["part"]
        if   part[0:5]=="kwire":
            return KwireComponents.KwireComponent.factory(ref,comp)
        elif part=="device:LED_Small": 
            return  DeviceLED_Small(ref,comp)
        elif part=="conn:CONN_01X02": 
            return  ConnCONN_01x02(ref,comp)
        elif part[0:4]=="74xx"     : 
            return SN74xx.SN74xx_4x2gate.factory(ref,comp)
        else:
            print("No model found for "+part)
            return Component(ref,comp)
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
        print("No model defined for component "+self.ref+" ("+self.part+")")
        pass
    def __str__(self):
        return self.__dict__.__str__()
    
class DeviceLED_Small(Component):
    def __str__(self):
        return self.__dict__.__str__()
    def __init__(self,ref,comp):
        self.A=False
        super().__init__(ref,comp)
    def setInput(self,context,time,pin,value):
        if value==self.A: return
        self.A=value
        print(self.ref,self.A)
        
class ConnCONN_01x02(Component):
    def __str__(self):
        return self.__dict__.__str__()
    def __init__(self,ref,comp):
        self.Y=[None]*(8+1)
        super().__init__(ref,comp)
        self.addOutput("1","GND")
        self.addOutput("2","VCC")
    def setInput(self,context,time,pin,value):
        if value==self.Y[int(pin)]: return
        self.Y[int(pin)]=value
        context.addEvent(time+self.propDelay,self.outputs[pin],value)
        print(self.ref,self.Y)
        
class Event:
    def __init__(self,time,net,value):
        self.time=time
        self.net=net
        self.value=value
    def __lt__(self,other):
        return self.time<other.time
    def __str__(self):
        return self.__dict__.__str__()

class Context:
    def __init__(self,infn):
        inf = antlr4.FileStream(infn)
        lexer = NetlistFileLexer(inf)
        stream = antlr4.CommonTokenStream(lexer)
        parser = NetlistFileParser(stream)
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
                if pin["type"]=="output" or pin["type"]=="3state":
                    comp.addOutput(node["pin"],name)
                    net["nodes"].remove(node)
                    
    def dispatchEvent(self):
        if len(self.eventList)==0:
            return False
        #Find the event with the earliest time
        event=heapq.heappop(self.eventList)
        print("Executing event: ",event)
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
    def __str__(self):
        return self.__dict__.__str__()
    
def toggleGrayCode(pins,tfirst):
    code=0
    result=[]
    for i in range(1,2**len(pins)+1):
        for j in range(32):
            if i and 2<<j:
                minbit=j
                break
        oldbit=(code>>j) and 1
        newbit=not oldbit
        code=(code and not (1<<minbit)) or (newbit<<minbit)
        result.append({"time":tfirst+i,"ref":pins[newbit]["ref"],"pin":pins[newbit]["pin"],"value":newbit==1})
    return result

if __name__ == '__main__':
    context=Context('kicad/xordriver.net')
    inputs=[{"time": 0,"ref":"J101","pin": "1","value":False}, #Turn on ground
            {"time": 0,"ref":"J101","pin": "2","value":True }]  #Turn on power
    inputs=inputs+toggleGrayCode([
            {"ref":"S101","pin":"16"},
            {"ref":"S102","pin":"16"}],0)
    for inp in inputs:
        print(inp)
    print(" T|B A|Y")
    print("--+---+-")
    for inp in inputs:
        context.components[inp["ref"]].setInput(context,inp["time"],inp["pin"],inp["value"])
        while context.dispatchEvent():
            print("Remaining events: ")
            for event in context.eventList:
                print(event)
            pass
        print("%2d|%s %s|%s" % (context.time,
              str(context.components["S102"].Y[1])[0],
              str(context.components["S101"].Y[1])[0],
              str(context.components["D201"].A   )[0]))
