##########################Working code for a 20-level V-cycle###################

from dolfin import *
import petsc4py.PETSc as pc
import re

##define the number of timesteps
L = 1
W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

source_concentration = Constant(.790)

##define the mesh and function space
#mesh = Mesh('GAMer_mesh_smoothed.xml')
mesh = BoxMesh(Point(0,0,0), Point(L,W,W), 10, 3, 3)#('/Users/mvhsan/Papers/mesh_dendritic_spine/final/tt-sr-mit.tet_mesh_smoothed_ODT.xml')#UnitSquareMesh(8, 8)
V = VectorFunctionSpace(mesh, 'P', 1)

##define the boundary conditions (here set to 1 + x^2 + 2y^2)
tol = 1E-14
def clamped_boundary(x, on_boundary):
	return on_boundary and x[0] < tol

bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)

#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
#bc = DirichletBC(V, u_D, 'on_boundary')

##define functions and differential equation
def epsilon(u):
	return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u):
	return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0, 0, -rho*g))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds
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
#x.create(pc.COMM_WORLD)
#x.setSizes(V.dim())
#x.setUp()
#x.assemble()

#if n == 0:
x_1 = pc.Vec()
x_1.create(pc.COMM_WORLD)
x_1.setSizes(V.dim())
x_1.setUp()

x_1.zeroEntries()
x_1.assemble()
#else:
#	x_1 = x1_set[0]

#P = decimation(mesh)
# P_set = {}
# for p in range(3):
# 	with open('/Users/mvhsan/Papers/fetk/mc/examples/simpleMCsh/pro_' + str(p) + '.m') as file:
# 		for line in file:
# 			if 'dim=' in line:
# 				dim1 = re.findall('([0-9]+)x[0-9]+', line)[0]
# 				dim2 = re.findall('[0-9]+x([0-9]+)', line)[0]
# 				break
# 		P = pc.Mat()
# 		P.create(pc.COMM_WORLD)
# 		P.setSizes([int(dim1), int(dim2)])
# 		P.setUp()
# 		for line in file:
# 			if 'e+' or 'e-' in line:
# 				i = re.findall('  ([0-9]+)  ', line)[0]
# 				j = re.findall('  ([0-9]+)  ', line)[1]
# 				val = re.findall('((?:\-|)[0-9]+\.[0-9]+)', line)[0]
# 				power = re.findall('e(?:\-|)([0-9]+)', line)[0]
# 				P[i, j] = val * pow(10, power)
# 	P.assemble()
# 	P_set[2-p] = P
P_set = {}
P = pc.Mat()
P.create(pc.COMM_WORLD)
P.setSizes(528)
P.setUp()
P.assemble()
y = []
for n in range(528):
    y.append(n)
P.zeroRowsColumns(y)
for i in range(3):
	P_set[i] = P


tol = 1e-12

x_set = {}
A_set = {}
b_set = {}
x1_set = {}
r_set = {}
e_set = {}

A_set[0] = A
b_set[0] = b
x1_set[0] = x_1

while (b - A*(x1_set[0])).norm() > tol:
	print((b-A*x1_set[0]).norm())
	for i in range(0, 3):
		x_set[i] = pc.Vec()
		if i == 0:
			(x1_set[i]).copy(result=(x_set[i]))
		
		##create discretization w/algebraic multigrid
		P_T = pc.Mat()
		P_set[i].transpose(out=P_T) #remember P_set[i] in real life
		
		if i != 0:
			transitional = P*A_set[i-1]
			A_set[i] = transitional*P_T
			#remember P_set[i] in real life
			#b_set[i] = P_T*b[i-1]	
			x_set[i].create(pc.COMM_WORLD)
			x_set[i].setSizes(A_set[i].size[1])
			x_set[i].setUp()

		##do Gauss-Seidel
		A_set[i].SOR(b_set[i], x_set[i], its=50, lits=100)
		r_set[i] = b_set[i] - A_set[i]*x_set[i]
		#print('r size ', str(r_set[i].size) + '\nP_T size ', str(P_T.size))
		b_set[i+1] = P_T*r_set[i]

	A_set[3] = P_set[2]*A_set[2] #remember P_set[18]
	
	p = pc.PC()
	p.create(pc.COMM_WORLD)
	p.setType(pc.PC.Type.LU)
	p.setOperators(A=A_set[3])
	p.setUp()

	k = pc.KSP()
	k.create(pc.COMM_WORLD)
	k.setPC(p)

	e_set[3] = pc.Vec()
	e_set[3].create(pc.COMM_WORLD)
	e_set[3].setSizes(A_set[3].size[1])
	e_set[3].setUp()
	
	k.solve(b_set[3], e_set[3])

	for i in range(3, 0, -1):
		e_set[i-1] = P*e_set[i]
		x1_set[i-1] = e_set[i-1] + x_set[i-1]
		A_set[i-1].SOR(b_set[i-1], x1_set[i-1], its=50, lits=100)
	
	print((b_set[0]-A_set[0]*x1_set[0]).norm())
	#soln[n] = x1_set[0]
#print(x1_set[0].getValues(y))
soln = PETScVector(x1_set[0])
outfile = File('testout.xml')
outfile << soln
meshout = File('mesh.pvd')
meshout << mesh

newmap = dof_to_vertex_map(V)

xmlValues = {}
mappedValues = {}

##putting the values in the output xml file into a dictionary mapping degree of freedom index to vector value
for line in open('testout.xml'):
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
with open('output.pvd', 'w+') as output:
    for line in open('mesh.pvd'):
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
