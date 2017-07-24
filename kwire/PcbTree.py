import antlr4
import heapq
from PcbFile.PcbFileLexer   import PcbFileLexer 
from PcbFile.PcbFileParser  import PcbFileParser
from PcbFile.PcbFileVisitor import PcbFileVisitor
 
class PcbTree(PcbFileVisitor):
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
    def visitString(self, ctx:PcbFileParser.StringContext):
        result=ctx.children[0].symbol.text
        if result[0]=='"':
            result=result[1:-1]
        return result
    def visitBoringNode(self, ctx:PcbFileParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:result[1]}
    def visitBoringIntNode(self, ctx:PcbFileParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:int(result[1])}
    def visitOptionalBoringNode(self, ctx:PcbFileParser.StringContext):
        result=self.visitChildren(ctx)
        if len(result) < 3:
            return {result[0][1:]:""}
        return {result[0][1:]:result[1]}
    #Nontrivial Nonterminals
    def visitPcb(self,ctx:PcbFileParser.StringContext):
        result=self.visitChildren(ctx)
        return result[0]
    def visitKicad_pcb(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitComponents(self,ctx:PcbFileParser.StringContext):
        children=self.visitChildren(ctx)
        components={}
        for component in children[1:-1]:
            ref=component["ref"]
            del component["ref"]
            components[ref]=component
        return {"components":components}
    def visitComp(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitLibsource(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libsource":dicts,"part":dicts["lib"]+":"+dicts["part"]}
    def visitLibparts(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libparts":dicts}
    def visitLibpart(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {dicts["lib"]+":"+dicts["part"]:dicts}
    def visitPins(self,ctx:PcbFileParser.StringContext):
        children=self.visitChildren(ctx)
        pins={}
        for pin in children[1:-1]:
            num=pin["num"]
            del pin["num"]
            pins[num]=pin
        return {"pins":pins}
    def visitPin(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitFields(self,ctx:PcbFileParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"fields":dicts}
    def visitField(self,ctx:PcbFileParser.StringContext):
        (dicts,nondicts)=self.sep(self.visitChildren(ctx))
        return {dicts["name"]:nondicts[1]}
    def visitNets(self,ctx:PcbFileParser.StringContext):
        children=self.visitChildren(ctx)
        nets={}
        for net in children[1:-1]:
            name=net["name"]
            del net["name"]
            nets[name]=net
        return {"nets":nets}
    def visitNet(self,ctx:PcbFileParser.StringContext):
        (dicts,nondicts)=self.sep(self.visitChildren(ctx))
        nodes=[]
        for node in nondicts[1:-1]:
            (nodedict,_)=self.sep(node)
            nodes.append(nodedict)
        dicts.update({"nodes":nodes})
        return dicts
    
    #Don't care nonterminals
    def visitDesign(self,ctx:PcbFileParser.StringContext): return None
    #Boring nodes
    def visitAlias      (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitCode       (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitDescription(self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitDocs       (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitFootprint  (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitFp         (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitLib        (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitLogical    (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitName       (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitNames      (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitNum        (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitNumber     (self,ctx:PcbFileParser.StringContext):return self.visitBoringIntNode(ctx)
    def visitPart       (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitPin__innode(self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitRef        (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitSource     (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitTool       (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamp     (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamps    (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitType_      (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitUri        (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitValue      (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    def visitVersion    (self,ctx:PcbFileParser.StringContext):return self.visitBoringNode(ctx)
    #Optional boring nodes
    def visitCompany(self,ctx:PcbFileParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitDate   (self,ctx:PcbFileParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitRev    (self,ctx:PcbFileParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitTitle  (self,ctx:PcbFileParser.StringContext):return self.visitOptionalBoringNode(ctx)

