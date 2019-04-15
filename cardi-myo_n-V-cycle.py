##########################Working code for a 20-level V-cycle###################

from dolfin import *
import petsc4py.PETSc as pc
import re

##where all output files go
filepref = 'cardi-myo_solution/'

##define number of multigrid levels to do
num_levels = 3

##define the number of timesteps
T = 1000
num_steps = 100 ##time-dependent
dt = T/num_steps
D = 79
k = Constant(dt)
D = Constant(D)

source_concentration = Constant(.790)

##define the mesh and function space
#mesh = Mesh('GAMer_mesh_smoothed.xml')
mesh = Mesh('/Users/mvhsan/Papers/mesh_dendritic_spine/final/tt-sr-mit.tet_mesh_smoothed_ODT.xml')#UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

meshxml = File(filepref + 'mesh.xml')
meshxml << mesh

Simplices = {}
#with open('tetmeshFromBlender_uninverted.xml') as file:
with open(filepref + 'mesh.xml') as file:
    for line in file:
        if '<tetrahedron ' in line:
            a_ = re.findall('index="([0-9]*)', line)
            a = a_[0]  ##or whatever the syntax is
            v0_ = re.findall('v0="([0-9]*)', line)	##need to fix syntax here
            v1_ = re.findall('v1="([0-9]*)', line)
            v2_ = re.findall('v2="([0-9]*)', line)
            v3_ = re.findall('v3="([0-9]*)', line)
            v0 = v0_[0]
            v1 = v1_[0]
            v2 = v2_[0]
            v3 = v3_[0]
            Simplices[a] = [v0, v1, v2, v3]
        if '<triangle ' in line:
            a_ = re.findall('index="([0-9]*)', line)
            a = a_[0]  ##or whatever the syntax is
            v0_ = re.findall('v0="([0-9]*)', line)	##need to fix syntax here
            v1_ = re.findall('v1="([0-9]*)', line)
            v2_ = re.findall('v2="([0-9]*)', line)
            v0 = v0_[0]
            v1 = v1_[0]
            v2 = v2_[0]
            Simplices[a] = [v0, v1, v2]

newmap = dof_to_vertex_map(V)

##define the boundary conditions (here set to 1 + x^2 + 2y^2)
facet_function = MeshFunction('size_t', mesh, 2, mesh.domains())

#bc_srrelease = DirichletBC(V, source_concentration, facet_function, 21)
bc = DirichletBC(V, source_concentration, facet_function, 21)
#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
#bc = DirichletBC(V, u_D, 'on_boundary')

##define functions and differential equation
phi = TrialFunction(V)
phi_n = Function(V)
v = TestFunction(V)
f = Constant(0)
sol = Function(V)
F = phi*v*dx + k*D*dot(grad(phi), grad(v))*dx - (phi_n + k*f)*v*dx
a, L = lhs(F), rhs(F)
#u = TrialFunction(V)
#v = TestFunction(V)
#f = Constant(-6.0)
#a = dot(grad(u), grad(v))*dx
#L = f*v*dx

solns = {}

num_P = num_levels

P = pc.Mat()
P.create(pc.COMM_WORLD)
P.setSizes(len(newmap))
P.setUp()
P.assemble()
y = []
for n in range(len(newmap)):
    y.append(n)
P.zeroRowsColumns(y)

P_set = {}
for p in range(num_P):
	P_set[p] = P

#P = decimation(mesh)
# P_set = {}
# num_P = num_levels - 1
# for p in range(num_P):
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
# 	P_set[num_P - 1 - p] = P

tol = 1e-8

t = 0
for n in range(num_steps):
	t += k
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

	if n == 0:
		x_1 = pc.Vec()
		x_1.create(pc.COMM_WORLD)
		x_1.setSizes(V.dim())
		x_1.setUp()

		x_1.zeroEntries()
		x_1.assemble()
	else:
		x_1 = x1_set[0]

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
		for i in range(num_levels):
			x_set[i] = pc.Vec()
			if i == 0:
				(x1_set[i]).copy(result=(x_set[i]))
		
			##create discretization w/algebraic multigrid
			P_T = pc.Mat()
			P_set[i].transpose(out=P_T) #remember P_set[i] in real life
		
			if i != 0:
				transitional = P_set[i-1]*A_set[i-1]
				A_set[i] = transitional*P_T
				#remember P_set[i] in real life
				#b_set[i] = P_T*b[i-1]	
				x_set[i].create(pc.COMM_WORLD)
				x_set[i].setSizes(A_set[i].size[1])
				x_set[i].setUp()

			##do Gauss-Seidel
			A_set[i].SOR(b_set[i], x_set[i], its=10, lits=10)
			r_set[i] = b_set[i] - A_set[i]*x_set[i]
			b_set[i+1] = P_T*r_set[i]

		#P_T = pc.Mat()
		#P_set[num_levels - 1].transpose(out=P_T)
		transitional = P_set[num_levels - 1]*A_set[num_levels - 1] #remember P_set[18]
		A_set[num_levels] = transitional*P_T
		
		p = pc.PC()
		p.create(pc.COMM_WORLD)
		p.setType(pc.PC.Type.HYPRE)
		#p.setType(pc.PC.Type.LU)
		p.setOperators(A=A_set[num_levels])
		p.setUp()

		ksp = pc.KSP()
		ksp.create(pc.COMM_WORLD)
		ksp.setPC(p)
		ksp.setType(pc.KSP.Type.CG)

		e_set[num_levels] = pc.Vec()
		e_set[num_levels].create(pc.COMM_WORLD)
		e_set[num_levels].setSizes(A_set[num_levels].size[1])
		e_set[num_levels].setUp()
		
		ksp.solve(b_set[num_levels], e_set[num_levels])

		for i in range(num_levels, 0, -1):
			e_set[i-1] = P_set[i - 1]*e_set[i]
			x1_set[i-1] = e_set[i-1] + x_set[i-1]
			A_set[i-1].SOR(b_set[i-1], x1_set[i-1], its=10, lits=10)
		
		print((b_set[0]-A_set[0]*x1_set[0]).norm())
	solns[n] = PETScVector(x1_set[0])
	sol = x1_set[0]

	outfile = File(filepref + 'function' + str(n) + '.xml')
	outfile << solns[n]

	Cell_dofs = {}
	for simp in Simplices:
		for dof in range(len(Simplices[simp])):
			vert = Simplices[simp][dof]
			try:
				Cell_dofs[vert]
			except:
				Cell_dofs[vert] = [simp, dof]

	with open(filepref + 'intermediary.xml', 'w+') as file:
		file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
		for index in range(sol.size):
			file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[index])][1]) + '" />\n')
			#file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
		file.write('  </function_data>\n</dolfin>')

	step = Function(V, filepref + 'intermediary.xml')
	phi_n.assign(step)
	#print(x1_set[0].getValues(y))


#sol = x1_set[0]
#soln = PETScVector(x1_set[0])
#outfile = File('testout.xml')
#meshfile = File('circlemesh.xml')
#meshfile << mesh
#outfile << soln
meshout = File(filepref + 'mesh.pvd')
meshout << mesh

xmlValues = {}
mappedValues = {}

##putting the values in the output xml file into a dictionary mapping degree of freedom index to vector value
for n in range(num_steps):
	xmlValues[n] = {}
	for line in open(filepref + 'function' + str(n) + '.xml'):
	    if('value=' in line):
	        m1 = re.search('(?<=value=")(\-|)[0-9]*\.[0-9]*e(\+|\-)[0-9]*', line)
	        m2 = re.search('(?<=index=")[0-9]*', line)
	        string1 = m1.group(0)
	        string2 = m2.group(0)
	        xmlValues[n][string2] = string1

for n in range(num_steps):
	mappedValues[n] = {}
	for i in range(0, len(newmap)):
	    mappedValues[n][str(newmap[i])] = xmlValues[n][str(i)]

t = 0
n = 0
#num_steps = 1
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
    with open(filepref + 'output_poisson' + fileid + '.vtu', 'w+') as output:
        #for line in open('mesh' + fileid + '.vtu'):
        for line in open(filepref + 'mesh000000.vtu'):
            if('</Cells>' in line):
                line = line + '<PointData  Scalars="f_8">\n<DataArray  type="Float64"  Name="f_8"  format="ascii">'
                #line = line + '<PointData Vectors="f_9">\n<DataArray  type="Float64"  Name="f_9"  NumberOfComponents="3" format="ascii">'
                for i in range(0, len(newmap)):
                    line = line + mappedValues[n][str(i)]
                    if(i != (len(newmap) - 1)):
                        line = line + '  '
                    else:
                        line = line + '</DataArray>\n</PointData>\n'
            output.write(line)
with open(filepref + 'output_poisson.pvd', 'w+') as output:
    for line in open(filepref + 'mesh.pvd'):
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
