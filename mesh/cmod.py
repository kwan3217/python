from OpenGLContext import testingcontext
BaseContext = testingcontext.getInteractive()
from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGLContext.arrays import *
from OpenGL.GL import shaders
from CMOD2Mesh2 import CMOD2Mesh2
import numpy as np

class TestContext( BaseContext ):
    """This shader just passes gl_Color from an input array to
    the fragment shader, which interpolates the values across the
    face (via a "varying" data type).
    """
    def __init__(self,mesh):
        self.mesh=mesh #Have to do this first, as BaseContext.__init() calls OnInit()
        BaseContext.__init__(self)
    def OnInit( self ):
        """Initialize the context once we have a valid OpenGL environ"""
        vertex = shaders.compileShader(
            """
            varying vec4 vertex_color;
            void main() {
                gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
                vertex_color = gl_Color;
            }""",GL_VERTEX_SHADER)
        fragment = shaders.compileShader("""
            varying vec4 vertex_color;
            void main() {
                gl_FragColor = vertex_color;
            }""",GL_FRAGMENT_SHADER)
        self.shader = shaders.compileProgram(vertex,fragment)
        #Geometry
        self.n_triangles=0
        for m in self.mesh.meshes:
            self.n_triangles+=len(m.groups[0].data)//3
        self.vbo=np.zeros((self.n_triangles*3,6),'f')
        i=0
        for mesh,matl in zip(self.mesh.meshes,self.mesh.materials):
            for v in mesh.groups[0].data:
                self.vbo[i,0:3]=mesh.V.data[CMOD2Mesh2.VertexSemantic.Position][v]*0.01
                self.vbo[i,3:6]=(matl.diffuse.red,matl.diffuse.green,matl.diffuse.blue)
                i+=1
        self.vbo = vbo.VBO(self.vbo)
        self.frame_count=0
    def Render(self, mode):
        self.frame_count += 1
        """Render the geometry for the scene."""
        BaseContext.Render(self, mode)
        glUseProgram(self.shader)
        try:
            self.vbo.bind()
            try:
                glEnableClientState(GL_VERTEX_ARRAY)
                glEnableClientState(GL_COLOR_ARRAY)
                glVertexPointer(3, GL_FLOAT, 24, self.vbo)
                glColorPointer(3, GL_FLOAT, 24, self.vbo + 12)
                glDrawArrays(GL_TRIANGLES, 0, self.n_triangles*3)
            finally:
                self.vbo.unbind()
                glDisableClientState(GL_VERTEX_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
        finally:
            shaders.glUseProgram( 0 )

if __name__ == "__main__":
    TestContext.ContextMainLoop(CMOD2Mesh2(open("Data/Mesh/voyager.cmod", "rb")))