import antlr4
from NetlistLexer   import NetlistLexer
from NetlistParser  import NetlistParser
from NetlistVisitor import NetlistVisitor
 
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
        return {"components":children[1:-1]}
    def visitComp(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitLibsource(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return {"libsource":dicts}
    def visitLibparts(self,ctx:NetlistParser.StringContext):
        children=self.visitChildren(ctx)
        return {"libparts":children[1:-1]}
    def visitLibpart(self,ctx:NetlistParser.StringContext):
        (dicts,_)=self.sep(self.visitChildren(ctx))
        return dicts
    def visitPins(self,ctx:NetlistParser.StringContext):
        children=self.visitChildren(ctx)
        return {"pins":children[1:-1]}
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
        return {"nets":children[1:-1]}
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

def parseNetlist(infn):     
    inf = antlr4.FileStream(infn)
    lexer = NetlistLexer(inf)
    stream = antlr4.CommonTokenStream(lexer)
    parser = NetlistParser(stream)
    tree = parser.netlist()
    visitor=NetlistTree()
    return visitor.visit(tree)

if __name__ == '__main__':
    netlist=parseNetlist('../../kicad/kwire/8bitadder.net')
    for comp in netlist["components"]:
        print(comp)
    for net in netlist["nets"]:
        print(net)
