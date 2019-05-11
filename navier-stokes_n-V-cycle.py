##########################Working code for a 20-level V-cycle###################

from dolfin import *
from mshr import *
import numpy as np
import petsc4py.PETSc as pc
import re

filepref = 'navier-stokes_solution/'
num_levels = 3

##define the number of timesteps
T = 5.0            # final time
num_steps = 5000   # number of time steps
dt = T / num_steps # time step size
mu = 0.001         # dynamic viscosity
rho = 1            # density

# Create mesh
channel = Rectangle(Point(0, 0), Point(2.2, 0.41))
cylinder = Circle(Point(0.2, 0.2), 0.05)
domain = channel - cylinder
mesh = generate_mesh(domain, 64)

# Define function spaces
V = VectorFunctionSpace(mesh, 'P', 1)
Q = FunctionSpace(mesh, 'P', 1)

# Define boundaries
inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 2.2)'
walls    = 'near(x[1], 0) || near(x[1], 0.41)'
cylinder = 'on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3'

# Define inflow profile
inflow_profile = ('4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0')

# Define boundary conditions
bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), inflow)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]


##define the mesh and function space
#mesh = Mesh('GAMer_mesh_smoothed.xml')

meshxml = File(filepref + 'mesh_ns.xml')
meshxml << mesh

Simplices = {}
#with open('tetmeshFromBlender_uninverted.xml') as file:
with open(filepref + 'mesh_ns.xml') as file:
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

Cell_dofs = {}
for simp in Simplices:
    for dof in range(len(Simplices[simp])):
        vert = Simplices[simp][dof]
        try:
            Cell_dofs[vert]
        except:
            Cell_dofs[vert] = [simp, dof]

##define the boundary conditions (here set to 1 + x^2 + 2y^2)


#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
#bc = DirichletBC(V, u_D, 'on_boundary')

##define functions and differential equation
# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n_  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx \
   + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
   + inner(sigma(U, p_n), epsilon(v))*dx \
   + dot(p_n*n_, v)*ds - dot(mu*nabla_grad(U)*n_, v)*ds \
   - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx


xdmffile_u = File(filepref + 'velocity.pvd')#XDMFFile(filepref + 'navier_stokes_cylinder/velocity.xdmf')
xdmffile_p = File(filepref + 'pressure.pvd')#XDMFFile(filepref + 'navier_stokes_cylinder/pressure.xdmf')


newmap = dof_to_vertex_map(V)
newmapQ = dof_to_vertex_map(Q)

solns = {}

num_P = num_levels

P = pc.Mat()
P.create(pc.COMM_WORLD)
#P.setSizes(len(newmapV))
P.setSizes(V.dim())
P.setUp()
P.assemble()
y = []
#for n in range(len(newmapV)):
for n in range(V.dim()):
    y.append(n)
P.zeroRowsColumns(y)

P_set = {}
for p_itr in range(num_P):
    P_set[p_itr] = P

P = pc.Mat()
P.create(pc.COMM_WORLD)
#P.setSizes(len(newmapQ))
P.setSizes(Q.dim())
P.setUp()
P.assemble()
y = []
#for n in range(len(newmapQ)):
for n in range(Q.dim()):
    y.append(n)
P.zeroRowsColumns(y)

P_setQ = {}
for p_itr in range(num_P):
    P_setQ[p_itr] = P

tol = 1e-10

A_petsc = PETScMatrix()
assemble(a1, tensor=A_petsc)
[bc.apply(A_petsc) for bc in bcu]
A1 = A_petsc.mat()

A_petsc = PETScMatrix()
assemble(a2, tensor=A_petsc)
[bc.apply(A_petsc) for bc in bcp]
A2 = A_petsc.mat()

A_petsc = PETScMatrix()
assemble(a3, tensor=A_petsc)
[bc.apply(A_petsc) for bc in bcu]
A3 = A_petsc.mat()


t = 0
for n in range(num_steps):
    t += dt

    print('starting timestep ' + str(n) + ' out of ' + str(num_steps))

    ##discretization into linear algebra matrix equation Ax=b
    A1_petsc = PETScMatrix()
    assemble(a1, tensor=A1_petsc)
    b_petsc = PETScVector()
    assemble(L1, tensor=b_petsc)
    [bc.apply(A1_petsc, b_petsc) for bc in bcu]
    #[bc.apply(b_petsc) for bc in bcu]

    ##convert to petsc4py objects
    A = A1_petsc.mat()
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
            A_set[i].SOR(b_set[i], x_set[i], its=20, lits=20)
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
            A_set[i-1].SOR(b_set[i-1], x1_set[i-1], its=20, lits=20)
        
        print((b_set[0]-A_set[0]*x1_set[0]).norm())
    solns[n] = PETScVector(x1_set[0])
    sol = x1_set[0]

    # Cell_dofs = {}
    # for simp in Simplices:
    #     for dof in range(len(Simplices[simp])):
    #         vert = Simplices[simp][dof]
    #         try:
    #             Cell_dofs[vert]
    #         except:
    #             Cell_dofs[vert] = [simp, dof]

    with open(filepref + 'intermediary.xml', 'w+') as file:
        file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
        # for index in range(sol.size):
        #     #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/3-(newmap[index]%3)/3))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/3-(newmap[index]%3)/3))][1]) + '" />\n')
        #     file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/2-(newmap[index]%2)/2))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/2-(newmap[index]%2)/2))][1]) + '" />\n')
        #     #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
        count = 0
        for index in range(sol.size):
            #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/3-(newmap[index]%3)/3))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/3-(newmap[index]%3)/3))][1]) + '" />\n')
            if (count % 2) == 0:
                file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/2))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/2))][1]) + '" />\n')
            else:
                file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/2))][0]) + '" cell_dof_index="' + str(3 + Cell_dofs[str(int(newmap[index]/2))][1]) + '" />\n')
            #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
            count += 1
        file.write('  </function_data>\n</dolfin>')

    step = Function(V, filepref + 'intermediary.xml')
    u_.assign(step)

    # outfile = File(filepref + 'function' + str(n) + '.xml')
    # outfile << solns[n]

    # output = File(filepref + 'solution.pvd')
    # output << solns[n]

    # Cell_dofs = {}
    # for simp in Simplices:
    #     for dof in range(len(Simplices[simp])):
    #         vert = Simplices[simp][dof]
    #         try:
    #             Cell_dofs[vert]
    #         except:
    #             Cell_dofs[vert] = [simp, dof]

    # with open(filepref + 'intermediary.xml', 'w+') as file:
    #     file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
    #     for index in range(sol.size):
    #         file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[index])][1]) + '" />\n')
    #         #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
    #     file.write('  </function_data>\n</dolfin>')

    # step = Function(V, filepref + 'intermediary.xml')
    # phi_n.assign(step)

    A2_petsc = PETScMatrix()
    assemble(a2, tensor=A2_petsc)
    b_petsc = PETScVector()
    assemble(L2, tensor=b_petsc)
    [bc.apply(A2_petsc, b_petsc) for bc in bcp]
    #[bc.apply(b_petsc) for bc in bcp]


    ##convert to petsc4py objects
    A = A2_petsc.mat()
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
        x_1.setSizes(Q.dim())
        x_1.setUp()

        x_1.zeroEntries()
        x_1.assemble()
    else:
        x_1 = x1_setQ[0]

    x_set = {}
    A_set = {}
    b_set = {}
    x1_setQ = {}
    r_set = {}
    e_set = {}

    A_set[0] = A
    b_set[0] = b
    x1_setQ[0] = x_1

    while (b - A*(x1_setQ[0])).norm() > tol:
        print((b-A*x1_setQ[0]).norm())
        for i in range(num_levels):
            x_set[i] = pc.Vec()
            if i == 0:
                (x1_setQ[i]).copy(result=(x_set[i]))
        
            ##create discretization w/algebraic multigrid
            P_T = pc.Mat()
            P_setQ[i].transpose(out=P_T) #remember P_set[i] in real life
        
            if i != 0:
                transitional = P_setQ[i-1]*A_set[i-1]
                A_set[i] = transitional*P_T
                #remember P_set[i] in real life
                #b_set[i] = P_T*b[i-1]  
                x_set[i].create(pc.COMM_WORLD)
                x_set[i].setSizes(A_set[i].size[1])
                x_set[i].setUp()

            ##do Gauss-Seidel
            A_set[i].SOR(b_set[i], x_set[i], its=20, lits=20)
            r_set[i] = b_set[i] - A_set[i]*x_set[i]
            b_set[i+1] = P_T*r_set[i]

        #P_T = pc.Mat()
        #P_set[num_levels - 1].transpose(out=P_T)
        transitional = P_setQ[num_levels - 1]*A_set[num_levels - 1] #remember P_set[18]
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
            e_set[i-1] = P_setQ[i - 1]*e_set[i]
            x1_setQ[i-1] = e_set[i-1] + x_set[i-1]
            A_set[i-1].SOR(b_set[i-1], x1_setQ[i-1], its=20, lits=20)
        
        print((b_set[0]-A_set[0]*x1_setQ[0]).norm())
    solns[n] = PETScVector(x1_setQ[0])
    sol = x1_setQ[0]

    # Cell_dofs = {}
    # for simp in Simplices:
    #     for dof in range(len(Simplices[simp])):
    #         vert = Simplices[simp][dof]
    #         try:
    #             Cell_dofs[vert]
    #         except:
    #             Cell_dofs[vert] = [simp, dof]

    with open(filepref + 'intermediaryQ.xml', 'w+') as file:
        file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
        for index in range(sol.size):
            file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmapQ[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmapQ[index])][1]) + '" />\n')
            #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmapQ[index])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmapQ[index])][1]) + '" />\n')
            #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
        file.write('  </function_data>\n</dolfin>')


    step = Function(Q, filepref + 'intermediaryQ.xml')
    p_.assign(step)

    A3_petsc = PETScMatrix()
    assemble(a3, tensor=A3_petsc)
    b_petsc = PETScVector()
    assemble(L3, tensor=b_petsc)
    [bc.apply(A3_petsc, b_petsc) for bc in bcu]

    ##convert to petsc4py objects
    A = A3_petsc.mat()
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
            A_set[i].SOR(b_set[i], x_set[i], its=20, lits=20)
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
            A_set[i-1].SOR(b_set[i-1], x1_set[i-1], its=20, lits=20)
        
        print((b_set[0]-A_set[0]*x1_set[0]).norm())
    solns[n] = PETScVector(x1_set[0])
    sol = x1_set[0]

    # Cell_dofs = {}
    # for simp in Simplices:
    #     for dof in range(len(Simplices[simp])):
    #         vert = Simplices[simp][dof]
    #         try:
    #             Cell_dofs[vert]
    #         except:
    #             Cell_dofs[vert] = [simp, dof]

    with open(filepref + 'intermediary.xml', 'w+') as file:
        file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <function_data size="' + str(sol.size) + '">\n')
        # for index in range(sol.size):
        #     file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/2-(newmap[index]%2)/2))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/2-(newmap[index]%2)/2))][1]) + '" />\n')
        #     #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
        count = 0
        for index in range(sol.size):
            #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/3-(newmap[index]%3)/3))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/3-(newmap[index]%3)/3))][1]) + '" />\n')
            if (count % 2) == 0:
                file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/2))][0]) + '" cell_dof_index="' + str(Cell_dofs[str(int(newmap[index]/2))][1]) + '" />\n')
            else:
                file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(int(newmap[index]/2))][0]) + '" cell_dof_index="' + str(3 + Cell_dofs[str(int(newmap[index]/2))][1]) + '" />\n')
            #file.write('    <dof index="' + str(index) + '" value="' + str(sol[index]) + '" cell_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][0]) + '" cell_dof_index="' + str(Cell_dofs[str(newmap[(index/3-(index%3))])][1]) + '" />\n')
            count += 1
        file.write('  </function_data>\n</dolfin>')

    step = Function(V, filepref + 'intermediary.xml')
    u_.assign(step)

    xdmffile_u << u_, t#xdmffile_u.write(u_, t)
    xdmffile_p << p_, t#xdmffile_p.write(p_, t)

    u_n.assign(u_)
    p_n.assign(p_)

    #saved_x1 = x1_set[0]



##putting the values in the output xml file into a dictionary mapping degree of freedom index to vector value
# meshout = File(filepref + 'mesh.pvd')
# meshout << mesh

# xmlValues = {}
# mappedValues = {}

# ##putting the values in the output xml file into a dictionary mapping degree of freedom index to vector value
# for n in range(num_steps):
#     xmlValues[n] = {}
#     for line in open(filepref + 'function' + str(n) + '.xml'):
#         if('value=' in line):
#             m1 = re.search('(?<=value=")(\-|)[0-9]*\.[0-9]*e(\+|\-)[0-9]*', line)
#             m2 = re.search('(?<=index=")[0-9]*', line)
#             string1 = m1.group(0)
#             string2 = m2.group(0)
#             xmlValues[n][string2] = string1

# for n in range(num_steps):
#     mappedValues[n] = {}
#     for i in range(0, len(newmap)):
#         mappedValues[n][str(newmap[i])] = xmlValues[n][str(i)]
        
# t = 0
# n = 0
# num_steps = 1
# ##here we define a function to dynamically create output file names (output000000, output000001, etc.)
# def getfileid(n):
#     if(n < 10):
#         fileid = '00000' + str(n)
#     elif(n < 100):
#         fileid = '0000' + str(n)
#     elif(n < 1000):
#         fileid = '000' + str(n)
#     elif(n < 10000):
#         fileid = '00' + str(n)
#     elif(n < 100000):
#         fileid = '0' + str(n)
#     else:
#         fileid = str(n)
#     return fileid
# for n in range(num_steps):
#     fileid = getfileid(n)
#     with open(filepref + 'output_elastic' + fileid + '.vtu', 'w+') as output:
#         for line in open(filepref + 'mesh_elastic000000.vtu'):
#             if('</Cells>' in line):
#                 #line = line + '<PointData  Scalars="f_8">\n<DataArray  type="Float64"  Name="f_8"  format="ascii">'
#                 line = line + '<PointData Vectors="f_9">\n<DataArray  type="Float64"  Name="f_9"  NumberOfComponents="3" format="ascii">'
#                 for i in range(0, len(newmap)):
#                     line = line + mappedValues[str(i)]
#                     if(i != (len(newmap) - 1)):
#                         line = line + '  '
#                     else:
#                         line = line + '</DataArray>\n</PointData>\n'
#             output.write(line)
# with open(filepref + 'output_elastic.pvd', 'w+') as output:
#     for line in open(filepref + 'mesh_elastic.pvd'):
#         newline = ''
#         if('file=' in line):
#             lineArray = re.split('(file="[a-z]*)', line)
#             for n in range(0, num_steps):
#                 fileid = getfileid(n)
#                 m = re.search('(?<=[0-9]{6}).*', lineArray[2], flags=re.DOTALL)
#                 newline = lineArray[0] + 'file="output_elastic' + fileid + m.group(0)
#                 output.write(newline)
#         else:
#             newline = line
#             output.write(newline)
