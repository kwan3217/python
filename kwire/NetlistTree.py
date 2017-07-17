import sys
import antlr4
from NetlistLexer   import NetlistLexer
from NetlistParser  import NetlistParser
from NetlistVisitor import NetlistVisitor
 
class KeyPrinter(NetlistVisitor):
    def visitTerminal(self, node):
        return node.symbol.text
    def aggregateResult(self, aggregate, nextResult):
        if aggregate is None:
            return [nextResult]
        if isinstance(aggregate,dict):
            aggregate.update(nextResult)
        aggregate.append(nextResult)
        return aggregate
    def visitString(self, ctx:NetlistParser.StringContext):
        result=ctx.children[0].symbol.text
        if result[0]=='"':
            result=result[1:-2]
        return result
    def visitBoringNode(self, ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        return {result[0][1:]:result[1]}
    def visitOptionalBoringNode(self, ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        if len(result) < 3:
            return {result[0][1:]:""}
        return {result[0][1:]:result[1]}
    #Nontrivial Nonterminals
    def visitExport(self,ctx:NetlistParser.StringContext):
        result=self.visitChildren(ctx)
        return result
    def visitComp(self,ctx:NetlistParser.StringContext):
    #Don't care nonterminals
    def visitDesign(self,ctx:NetlistParser.StringContext): return None
    #Boring nodes
    def visitAlias(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitCode(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitDescription(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitDocs(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitFootprint(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitFp(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitLib(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitLogical(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitName(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitNames(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitNum(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitNumber(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitPart(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitPin__innode(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitRef(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitSource(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitTool(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamp(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitTstamps(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitType(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitUri(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitValue(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    def visitVersion(self,ctx:NetlistParser.StringContext):return self.visitBoringNode(ctx)
    #Optional boring nodes
    def visitCompany(self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitDate(self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitRev(self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)
    def visitTitle(self,ctx:NetlistParser.StringContext):return self.visitOptionalBoringNode(ctx)

     

def main(argv):
    input = antlr4.FileStream('../../kicad/Tomo/Tomo.net')
    lexer = NetlistLexer(input)
    stream = antlr4.CommonTokenStream(lexer)
    parser = NetlistParser(stream)
    tree = parser.netlist()
    visitor=KeyPrinter()
    result=visitor.visit(tree)
    print(result)

if __name__ == '__main__':
    main(sys.argv)
