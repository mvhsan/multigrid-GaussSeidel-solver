##########################Working code for a 20-level V-cycle###################

from dolfin import *
import petsc4py.PETSc as pc
import re

##define the number of timesteps

##define the mesh and function space
#mesh = Mesh('GAMer_mesh_smoothed.xml')
from mshr import *
bigSphere = Sphere(Point(0, 0, 0), 100)
sphere1 = Sphere(Point(0, 10, 1), 5)
sphere2 = Sphere(Point(20, 30, 20), 10)
domain = bigSphere - sphere1
domain = domain - sphere2
mesh = generate_mesh(domain, 64)
V = FunctionSpace(mesh, 'P', 1)

#scalar values for F(u)
a_s = 1
b_s = 1
c_s = 1
d_s = 1

filepref = 'einstein_constrained_solution/'

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

##define the boundary conditions (here set to 1 + x^2 + 2y^2)
#w_D = Expression('1 - pow(x[0], 2) - pow(x[1], 2)', degree = 1)
w_D = Constant(1.0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, w_D, boundary)

n_D = Constant(0.0)

newton_bc = DirichletBC(V, n_D, boundary)

## Function F(u)
def F_funct(u):
    return -div(grad(u)) + a_s*u + b_s*pow(u, 5) - c_s*pow(u, -7) - d_s*pow(u, -3)

#F_expr = Expression('pow(x, 3) + div(grad(x))', degree=1)

#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
#bc = DirichletBC(V, u_D, 'on_boundary')

##define functions and differential equation
u = TrialFunction(V)
v = TestFunction(V)
u_ = Function(V)
u_n = Function(V)
f = Constant(0)
F = dot(grad(u), grad(v))*dx + (a_s*u + b_s*pow(u, 5) - c_s*pow(u, -7) - d_s*pow(u, -3))*v*dx
F = action(F, u_)
L1 = f*v*dx
L2 = -F_funct(u_)
Jacob = derivative(F, u_, u)
a = dot(grad(u), grad(v))*dx
print(type(F))
#print(type(J))
#u = TrialFunction(V)
#v = TestFunction(V)
#f = Constant(-6.0)
#a = dot(grad(u), grad(v))*dx
#L = f*v*dx

##discretization into linear algebra matrix equation Ax=b
# J_petsc = PETScMatrix()
# assemble(Jacob, tensor=J_petsc)
# h_petsc = PETScVector()
# assemble(L2, tensor=h_petsc)
# bc.apply(J_petsc, h_petsc)

# ##convert to petsc4py objects
# J = J_petsc.mat()
# #A.assemble()
# h = h_petsc.vec()
# #b.assemble()

newmap = dof_to_vertex_map(V)
print(newmap)
print(len(newmap))

A_petsc = PETScMatrix()
assemble(a, tensor=A_petsc)
b_petsc = PETScVector()
assemble(L1, tensor=b_petsc)
bc.apply(A_petsc, b_petsc)

A = A_petsc.mat()
b = b_petsc.vec()

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
P.setSizes(b.size)
P.setUp()
P.assemble()
y = []
for n in range(b.size):
    y.append(n)
P.zeroRowsColumns(y)
for i in range(3):
	P_set[i] = P


tol = 0.001#0.0044#1e-12

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
		A_set[i].SOR(b_set[i], x_set[i], its=10, lits=10)
		r_set[i] = b_set[i] - A_set[i]*x_set[i]
		#print('r size ', str(r_set[i].size) + '\nP_T size ', str(P_T.size))
		b_set[i+1] = P_T*r_set[i]

	A_set[3] = P_set[2]*A_set[2] #remember P_set[18]
	
	p = pc.PC()
	p.create(pc.COMM_WORLD)
	p.setType(pc.PC.Type.NONE)
	p.setOperators(A=A_set[3])
	p.setUp()

	k = pc.KSP()
	k.create(pc.COMM_WORLD)
	k.setPC(p)
	k.setType(pc.KSP.Type.CG)

	e_set[3] = pc.Vec()
	e_set[3].create(pc.COMM_WORLD)
	e_set[3].setSizes(A_set[3].size[1])
	e_set[3].setUp()
	
	k.solve(b_set[3], e_set[3])

	for i in range(3, 0, -1):
		e_set[i-1] = P*e_set[i]
		x1_set[i-1] = e_set[i-1] + x_set[i-1]
		A_set[i-1].SOR(b_set[i-1], x1_set[i-1], its=10, lits=10)
	
	print((b_set[0]-A_set[0]*x1_set[0]).norm())
print((b_set[0]-A_set[0]*x1_set[0]).norm())

	#soln[n] = x1_set[0]
#print(x1_set[0].getValues(y))
sol = x1_set[0]
soln = PETScVector(x1_set[0])

Cell_dofs = {}
for simp in Simplices:
    for dof in range(len(Simplices[simp])):
        vert = Simplices[simp][dof]
        try:
            Cell_dofs[vert]
        except:
            Cell_dofs[vert] = [simp, dof]


with open(filepref + 'intermediary_1.xml', 'w+') as file:
    file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
    for index in range(sol.size):
        file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[index])][1]) + '" />\n')
    file.write('  </function_data>\n</dolfin>')

step = Function(V, filepref + 'intermediary_1.xml')
u_.assign(step)

testfile = File(filepref + 'testingoutput.pvd')
testfile << u_

#print(type(F_funct(u_)))
#print(F_funct(u_))
#print(dir(F_funct(u_)))
#print(F_funct(u_).evaluate())

#zero = div(grad(u_))
#try1 = project(zero, V)#interpolate(zero, V)
#something = Expression("-Z+pow(U, 3)", Z=zero, U=u_, degree=1)

#print(type(F_expr))
#print(F_expr)

#one2 = div(grad(u_))
#one3 = Expression('a', a=one2, degree=1)
#one1 = Expression('div(grad(U))', U=u_, degree=1)
#two2 = Expression('pow(U, 3)', U=u_, degree=1)
#print(one1+two2)

#one = div(grad(u_))
#def one(u):
#    return div(grad(u))
#two = pow(u_, 3)
#print(one)
#print(two)
#print(type(one))
#print(type(two))
#print(type(u_))
#three = interpolate(one(u_), V)
#huh = Expression('x[0] + 1', degree=1)

#one = interpolate(huh, V)

check = sqrt(inner(F_funct(u_), F_funct(u_)))
#check = Expression('sqrt(inner(-div(grad(u))+pow(u, 3), -div(grad(u))+pow(u,3)))', u=u_, degree=1)
#check1 = Expression('F', F=F_funct(u_), degree=1)
#check = Expression('sqrt(inner(F, F))', F=F_funct(u_), degree=1)
#check1 = Expression('c', c=check, degree=1)
print(check)
print(type(check))
#print(check._ufl_evaluate_scalar_())
#thing1 = project(check)
#print(check._ufl_evaluate_scalar_(u_))
print(type(V))
#thing = interpolate(check, V)#project(check, V)
thing = project(check, V, solver_type="cg", preconditioner_type="amg")
#print(type(thing))
print('magnitude of F(u): ', norm(thing))

while (norm(thing) > tol):

    J_petsc = PETScMatrix()
    Jacob = derivative(F, u_, u)#(F, u_, u)
    assemble(Jacob, tensor=J_petsc)
    #h_petsc = PETScVector()
    #trying = project(L2, V, solver_type="cg", preconditioner_type="amg")
    F_ = project(-F_funct(u_), V, solver_type="cg", preconditioner_type="amg")
    F_ = F_.vector()
    h_petsc = as_backend_type(F_)
    #F_ = F_.vec()
    print(type(h_petsc))
    #assemble(L2, tensor=h_petsc)
    newton_bc.apply(J_petsc, h_petsc)

    ##convert to petsc4py objects
    J = J_petsc.mat()
    #A.assemble()
    h = h_petsc.vec()

    x = pc.Vec()
    #x.create(pc.COMM_WORLD)
    #x.setSizes(V.dim())
    #x.setUp()
    #x.assemble()

    #if n == 0:
    x_1 = sol#pc.Vec()
    # x_1.create(pc.COMM_WORLD)
    # x_1.setSizes(V.dim())
    # x_1.setUp()

    # x_1.zeroEntries()
    # x_1.assemble()
    #else:
    #   x_1 = x1_set[0]

    #P = decimation(mesh)
    # P_set = {}
    # for p in range(3):
    #   with open('/Users/mvhsan/Papers/fetk/mc/examples/simpleMCsh/pro_' + str(p) + '.m') as file:
    #       for line in file:
    #           if 'dim=' in line:
    #               dim1 = re.findall('([0-9]+)x[0-9]+', line)[0]
    #               dim2 = re.findall('[0-9]+x([0-9]+)', line)[0]
    #               break
    #       P = pc.Mat()
    #       P.create(pc.COMM_WORLD)
    #       P.setSizes([int(dim1), int(dim2)])
    #       P.setUp()
    #       for line in file:
    #           if 'e+' or 'e-' in line:
    #               i = re.findall('  ([0-9]+)  ', line)[0]
    #               j = re.findall('  ([0-9]+)  ', line)[1]
    #               val = re.findall('((?:\-|)[0-9]+\.[0-9]+)', line)[0]
    #               power = re.findall('e(?:\-|)([0-9]+)', line)[0]
    #               P[i, j] = val * pow(10, power)
    #   P.assemble()
    #   P_set[2-p] = P
    P_set = {}
    P = pc.Mat()
    P.create(pc.COMM_WORLD)
    P.setSizes(b.size)
    P.setUp()
    P.assemble()
    y = []
    for n in range(b.size):
        y.append(n)
    P.zeroRowsColumns(y)
    for i in range(3):
        P_set[i] = P



    x_set = {}
    A_set = {}
    b_set = {}
    x1_set = {}
    r_set = {}
    e_set = {}

    A_set[0] = J
    b_set[0] = h
    x1_set[0] = x_1

    while (h - J*(x1_set[0])).norm() > tol:
        print((h-J*x1_set[0]).norm())
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
            A_set[i].SOR(b_set[i], x_set[i], its=10, lits=10)
            r_set[i] = b_set[i] - A_set[i]*x_set[i]
            #print('r size ', str(r_set[i].size) + '\nP_T size ', str(P_T.size))
            b_set[i+1] = P_T*r_set[i]

        A_set[3] = P_set[2]*A_set[2] #remember P_set[18]
        
        p = pc.PC()
        p.create(pc.COMM_WORLD)
        p.setType(pc.PC.Type.NONE)
        p.setOperators(A=A_set[3])
        p.setUp()

        k = pc.KSP()
        k.create(pc.COMM_WORLD)
        k.setPC(p)
        k.setType(pc.KSP.Type.CG)

        e_set[3] = pc.Vec()
        e_set[3].create(pc.COMM_WORLD)
        e_set[3].setSizes(A_set[3].size[1])
        e_set[3].setUp()
        
        k.solve(b_set[3], e_set[3])

        for i in range(3, 0, -1):
            e_set[i-1] = P*e_set[i]
            x1_set[i-1] = e_set[i-1] + x_set[i-1]
            A_set[i-1].SOR(b_set[i-1], x1_set[i-1], its=10, lits=10)
        
        print((b_set[0]-A_set[0]*x1_set[0]).norm())
    print((b_set[0]-A_set[0]*x1_set[0]).norm())

        #soln[n] = x1_set[0]
    #print(x1_set[0].getValues(y))
    sol = sol + x1_set[0]
    soln = PETScVector(x1_set[0])

    Cell_dofs = {}
    for simp in Simplices:
        for dof in range(len(Simplices[simp])):
            vert = Simplices[simp][dof]
            try:
                Cell_dofs[vert]
            except:
                Cell_dofs[vert] = [simp, dof]


    with open(filepref + 'intermediary_1.xml', 'w+') as file:
        file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
        for index in range(sol.size):
            file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[index])][1]) + '" />\n')
        file.write('  </function_data>\n</dolfin>')

    step = Function(V, filepref + 'intermediary_1.xml')
    u_.assign(step)

    check = sqrt(inner(F_funct(u_), F_funct(u_)))
    thing = project(check, V, solver_type="cg", preconditioner_type="amg")
    print('magnitude of F(u): ', norm(thing))






outfile = File(filepref + 'testout_einstein.xml')
#meshfile = File('circlemesh.xml')
#meshfile << mesh
outfile << soln
meshout = File(filepref + 'mesh_einstein.pvd')
meshout << mesh
#xdmffile = XDMFFile('testout.xdmf')
#xdmffile << soln

xmlValues = {}
mappedValues = {}

# Cell_dofs = {}
# for simp in Simplices:
# 	for dof in range(len(Simplices[simp])):
# 		vert = Simplices[simp][dof]
# 		try:
# 			Cell_dofs[vert]
# 		except:
# 			Cell_dofs[vert] = [simp, dof]


# with open(filepref + 'intermediary.xml', 'w+') as file:
# 	file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
# 	for index in range(sol.size):
# 		file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[index])][1]) + '" />\n')
# 	file.write('  </function_data>\n</dolfin>')


##putting the values in the output xml file into a dictionary mapping degree of freedom index to vector value
for line in open(filepref + 'testout_einstein.xml'):
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
    with open(filepref + 'output_einstein' + fileid + '.vtu', 'w+') as output:
        for line in open(filepref + 'mesh_einstein' + fileid + '.vtu'):
            if('</Cells>' in line):
                line = line + '<PointData  Scalars="f_8">\n<DataArray  type="Float64"  Name="f_8"  format="ascii">'
                for i in range(0, len(newmap)):
                    line = line + mappedValues[str(i)]
                    if(i != (len(newmap) - 1)):
                        line = line + '  '
                    else:
                        line = line + '</DataArray>\n</PointData>\n'
            output.write(line)
with open(filepref + 'output_einstein.pvd', 'w+') as output:
    for line in open(filepref + 'mesh_einstein.pvd'):
        newline = ''
        if('file=' in line):
            lineArray = re.split('(file="[a-z]*)', line)
            for n in range(0, num_steps):
                fileid = getfileid(n)
                m = re.search('(?<=[0-9]{6}).*', lineArray[2], flags=re.DOTALL)
                newline = lineArray[0] + 'file="output_einstein' + fileid + m.group(0)
                output.write(newline)
        else:
            newline = line
            output.write(newline)
