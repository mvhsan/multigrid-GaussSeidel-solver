##########################Working code for a 20-level V-cycle###################

from dolfin import *
import petsc4py.PETSc as pc
import re

##define the number of timesteps
beta = 8
R0 = 0.6
p = Expression('4*exp(-pow(beta, 2)*(pow(x[0], 2) + pow(x[1] - R0, 2)))', degree = 1, beta = beta, R0 = R0)

##define the mesh and function space
#mesh = Mesh('GAMer_mesh_smoothed.xml')
from mshr import *
domain = Circle(Point(0, 0), 1)
mesh = generate_mesh(domain, 64)
V = FunctionSpace(mesh, 'P', 1)

##define the boundary conditions (here set to 1 + x^2 + 2y^2)
w_D = Expression('1 - pow(x[0], 2) - pow(x[1], 2)', degree = 1)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, w_D, boundary)

#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
#bc = DirichletBC(V, u_D, 'on_boundary')

##define functions and differential equation
w = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(w), grad(v))*dx
L = p*v*dx
#u = TrialFunction(V)
#v = TestFunction(V)
#f = Constant(-6.0)
#a = dot(grad(u), grad(v))*dx
#L = f*v*dx

##discretization into linear algebra matrix equation Ax=b
A_petsc = PETScMatrix()
assemble(a, tensor=A_petsc)
b_petsc = PETScVector()
assemble(L, tensor=b_petsc)
bc.apply(A_petsc, b_petsc)

##convert to petsc4py objects
A = A_petsc.mat()
#A.assemble()
b = b_petsc.vec()
#b.assemble()

##initialize x, the output vector
x = pc.Vec()
x.create(pc.COMM_WORLD)
x.setSizes(b.size)
x.setUp()
#x.create(pc.COMM_WORLD)
#x.setSizes(V.dim())
#x.setUp()
#x.assemble()

o = pc.Options()
o.create()
o.setValue('-pc_hypre_boomeramg_tol', 1E-12)
o.setValue('-pc_gamg_coarse_eq_limit', 10000000)

PC = pc.PC()
PC.create(pc.COMM_WORLD)
PC.setOperators(A=A)
#PC.setType(pc.PC.Type.GAMG)#HYPRE)
PC.setType(pc.PC.Type.HYPRE)
#PC.setGAMGType('geometric')
#PC.setGAMGLevels(1)
PC.setUp()
CG = pc.KSP()
CG.create(pc.COMM_WORLD)
#CG.setType(pc.KSP.Type.RICHARDSON)#PREONLY)
CG.setType(pc.KSP.Type.RICHARDSON)
CG.setPC(PC)
CG.solve(b, x)

print(((b-A*x).norm()))

soln = PETScVector(x)
outfile = File('amg_testout.xml')
outfile << soln
meshout = File('amg_mesh.pvd')
meshout << mesh

newmap = dof_to_vertex_map(V)

xmlValues = {}
mappedValues = {}

##putting the values in the output xml file into a dictionary mapping degree of freedom index to vector value
for line in open('amg_testout.xml'):
    if('value=' in line):
        m1 = re.search('(?<=value=")(\-|)[0-9]*\.[0-9]*e(\+|\-)[0-9]*', line)
        m2 = re.search('(?<=index=")[0-9]*', line)
        string1 = m1.group(0)
        string2 = m2.group(0)
        xmlValues[string2] = string1

for i in range(0, len(newmap)):
    mappedValues[str(newmap[i])] = xmlValues[str(i)]

t = 0
n = 0
num_steps = 1
##here we define a function to dynamically create output file names (output000000, output000001, etc.)
def getfileid(n):
    if(n < 10):
        fileid = '00000' + str(n)
    elif(n < 100):
        fileid = '0000' + str(n)
    elif(n < 1000):
        fileid = '000' + str(n)
    elif(n < 10000):
        fileid = '00' + str(n)
    elif(n < 100000):
        fileid = '0' + str(n)
    else:
        fileid = str(n)
    return fileid
for n in range(num_steps):
    fileid = getfileid(n)
    with open('output' + fileid + '.vtu', 'w+') as output:
        for line in open('mesh' + fileid + '.vtu'):
            if('</Cells>' in line):
                line = line + '<PointData  Scalars="f_8">\n<DataArray  type="Float64"  Name="f_8"  format="ascii">'
                for i in range(0, len(newmap)):
                    line = line + mappedValues[str(i)]
                    if(i != (len(newmap) - 1)):
                        line = line + '  '
                    else:
                        line = line + '</DataArray>\n</PointData>\n'
            output.write(line)
with open('amg_output.pvd', 'w+') as output:
    for line in open('amg_mesh.pvd'):
        newline = ''
        if('file=' in line):
            lineArray = re.split('(file="[a-z]*)', line)
            for n in range(0, num_steps):
                fileid = getfileid(n)
                m = re.search('(?<=[0-9]{6}).*', lineArray[2], flags=re.DOTALL)
                newline = lineArray[0] + 'file="output' + fileid + m.group(0)
                output.write(newline)
        else:
            newline = line
            output.write(newline)
