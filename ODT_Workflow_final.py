#! /Users/mvhsan/anaconda3/envs/fenicsproject/bin/python
###Code for the ODT mesh smoothing workflow

from dolfin import *
import petsc4py.PETSc as pc
import re
import os
import time

Tetrahedrons_ = {}
#with open('tetmeshFromBlender_uninverted.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
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
            Tetrahedrons_[a] = [v0, v1, v2, v3]

Original_vertices = {}
#n = 0
N = 0
#with open('tetmeshFromBlender_uninverted.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
    for line in file:
        if '<vertices' in line:
            N = int((re.findall('size="([0-9]*)', line))[0])
            #print(N)
        if '<vertex' in line:
            index = (re.findall('index="([0-9]*)', line))[0]
            #x = (re.findall(' x="([0-9]*.[0-9]*|-[0-9]*.[0-9]*)', line))[0]
            #x = (re.findall(' x="((?:[0-9]*|\-[0-9]*).[0-9]+)', line))[0]
            x = (re.findall(' x="((?:\-|)(?:[0-9]*\.[0-9]*|[0-9]*))', line))[0]
            y = (re.findall(' y="((?:\-|)(?:[0-9]*\.[0-9]*|[0-9]*))', line))[0]
            z = (re.findall(' z="((?:\-|)(?:[0-9]*\.[0-9]*|[0-9]*))', line))[0]
            Original_vertices[index] = [x, y, z]

#print('n: ', n)
print('N: ', N)
print('original vertices length: ', len(Original_vertices))

#N = 24237
Nb = 0
P_dict = {}
P_dict[0] = {}
Intermediary = {}
Reverse_intermediary = {}
#with open('tetmeshFromBlender_uninverted.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
    num = 0
    itr = 0
    for line in file:
        if '<mesh_value_collection' in line:
            num = (re.findall('size="([0-9]*)', line))[0]
            #print(num)
            break
    n = 0
    itr = int(num)
    N_step = 1
    D = 0
    for line in file:
        if itr > 0:
            #D = D + 1
            if '<value' in line:
                D = D + 1
                #if 'value="1"' not in line:
                #if 'value="21"' in line:
                if 'value="1"' in line:
                    entity = (re.findall('local_entity="([0-9]*)', line))[0]
                    tet = (re.findall('cell_index="([0-9]*)', line))[0]
                    for vert in Tetrahedrons_[tet]:
                        if vert != entity:
                            if vert not in Intermediary:
                                P_dict[0][N - N_step] = Original_vertices[vert]
                                Intermediary[vert] = N-N_step
                                Reverse_intermediary[N-N_step] = vert
                                N_step = N_step + 1
                    #itr = itr - 1
                #itr = itr - 1
            if 'dim="3"' and 'mesh_value_collection' in line:
                break
    #print(len(Intermediary))
    N = N - N_step + 1
    Nb = N_step - 1
    itr = int(num)
    n = 0

for vert in Original_vertices:
    if vert not in Intermediary:
        P_dict[0][n] = Original_vertices[vert]
        Intermediary[vert] = n
        Reverse_intermediary[n] = vert
        n = n + 1
    
print('n: ', n)
#with open('/Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml') as file:
    #print(itr)
    #file.seek(0)
#    S = 0
#    for line in file:
#        if 'mesh_value_collection' in line:
#            print(line)
#            break
#    for line in file:
#        #S = S + 1
#        if itr > 0:
#            #S = S + 1
#            if '<value' in line:
#                S = S + 1
#                #if 'value="1"' in line: ##or whatever is interior marker
                #if 'value="30"' in line:
                #if 'value="1"' not in line:
                #    entity = (re.findall('local_entity="([0-9]*)', line))[0]
                #    tet = (re.findall('cell_index="([0-9]*)', line))[0]
                #    for vert in Tetrahedrons_[tet]:
                #        if vert != entity:
                #            if vert not in Intermediary:
                #                P_dict[0][n] = Original_vertices[vert]
                #                Intermediary[vert] = n
                #                Reverse_intermediary[n] = vert
                #                itr = itr - 1
                #                n = n + 1
                #entity = (re.findall('local_entity="([0-9]*)', line))[0]
                #tet = (re.findall('cell_index="([0-9]*)', line))[0]
                #for vert in Tetrahedrons_[tet]:
                #    if vert != entity:
                #        #S = S + 1
                #        if vert not in Intermediary:
                #            #print('here')
                #            P_dict[0][n] = Original_vertices[vert]
                #            Intermediary[vert] = n
                #            Reverse_intermediary[n] = vert
                #            #itr = itr - 1
                #            n = n + 1
                #itr = itr - 1
            #if 'dim="3"' and 'mesh_value_collection' in line:
            #    break
    #print(len(Intermediary))

#print('S: ', S)
print('D: ', D)
print(len(P_dict[0]))
with open('P_dict_test.txt', 'w+') as file:
    for i in P_dict[0]:
        line = str(i) + '\n'
        file.write(line)
with open('Intermediary_test.txt', 'w+') as file:
    for i in Intermediary:
        line = str(i) + '\n'
        file.write(line)
print('interior vertices: ', N)
print('boundary vertices: ', Nb)
print('total vertices: ', N+Nb)

Tetrahedrons = {}
#with open('tetmeshFromBlender_uninverted.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
    for line in file:
        if '<tetrahedron' in line:
            a_ = re.findall('index="([0-9]*)', line)
            a = int(a_[0])  ##or whatever the syntax is
            v0_ = re.findall('v0="([0-9]*)', line)	##need to fix syntax here
            v1_ = re.findall('v1="([0-9]*)', line)
            v2_ = re.findall('v2="([0-9]*)', line)
            v3_ = re.findall('v3="([0-9]*)', line)
            v0 = int(v0_[0])
            v1 = int(v1_[0])
            v2 = int(v2_[0])
            v3 = int(v3_[0])
            #Tetrahedrons[a] = [Reverse_intermediary[v0], Reverse_intermediary[v1], Reverse_intermediary[v2], Reverse_intermediary[v3]]
            Tetrahedrons[a] = [str(Intermediary[str(v0)]), str(Intermediary[str(v1)]), str(Intermediary[str(v2)]), str(Intermediary[str(v3)])]

print(len(Tetrahedrons))

Vertices = {}
#for vert in range(N):		## if I fix the vertex indexing from the creation of Tetrahedrons then I should be able to just say range(N)
#    Vertices[vert] = []
#    for tet in range(len(Tetrahedrons)):
#        if vert in Tetrahedrons[tet]:
#            Vertices[vert].append(tet)
for tet in range(len(Tetrahedrons)):
    for vert in Tetrahedrons[tet]:
        try:
            Vertices[vert]
        except:
            Vertices[vert] = []
        Vertices[vert].append(tet)

print('#####################################################')
print(Vertices['0'])
print(Intermediary['0'])
#print(Tetrahedrons_['78011'])
#print(Tetrahedrons_['78014'])
#print(Tetrahedrons_['78023'])
#print(Tetrahedrons_['78027'])
#print(Tetrahedrons_['78033'])
#print(Tetrahedrons_['78037'])
print('#####################################################')

def star(vert):	##takes in index of the vertex
    result = []
    for tet in Vertices[vert]:
        result.append(tet)
    return result

A_bar = pc.Mat() ##N is number of interior points, Nb that of boundary points
A_bar.create(pc.COMM_WORLD)
A_bar.setSizes([N+Nb, N+Nb])
A_bar.setUp()
#A = pc.Mat()
#A.create(pc.COMM_WORLD)
#A.setSizes([N, N])
#A.setUp()
#I = 0
global I
I = 0
rho = 1 ##or error density function
Eq = {}
h = {}

Tetrahedron_matrix = {}
t_n = 0
for tetrahedron in Tetrahedrons:
    Tet = pc.Mat()
    Tet.create(pc.COMM_WORLD)
    Tet.setSizes([4, 3])
    Tet.setUp()
    #Tet.assemble()
    m = 0
    for vertex in Tetrahedrons[tetrahedron]:
        n = 0
        for entry in P_dict[I][int(vertex)]:#[int(vertex)]:
            #print(P_dict[i][int(vertex)])
            #print(entry)
            Tet[m, n] = float(entry)
            n = n + 1
        m = m + 1
    Tetrahedron_matrix[t_n] = Tet
    t_n = t_n + 1
    
Tetrahedron_array = {}
t_n = 0
for tetrahedron in Tetrahedrons:
    Tet = {}
    m = 0
    for vertex in Tetrahedrons[tetrahedron]:
        Tet[m] = P_dict[I][int(vertex)]
        m = m + 1
    Tetrahedron_array[t_n] = Tet
    t_n = t_n + 1

Tetrahedron_matrix[0].assemble()
print(Tetrahedron_matrix[0].getValues([0, 1, 2, 3], [0, 1, 2]))
print(Tetrahedron_array[0])
print('########################################')

def volume(j):
    #p = pc.PC()
    #p.create(pc.COMM_WORLD)
    #p.setType(pc.PC.Type.LU)
    #p.setOperators(A=Tetrahedron_matrix[j])
    #p.setUp()
    
    #k = pc.KSP()
    #k.create(pc.COMM_WORLD)
    #k.setPC(p)
    
    VecA = Tetrahedron_array[j][0]
    VecB = Tetrahedron_array[j][1]
    VecC = Tetrahedron_array[j][2]
    VecD = Tetrahedron_array[j][3]
    
    A_diff = {}
    B_diff = {}
    C_diff = {}
    for n in range(3):
        A_diff[n] = float(VecA[n]) - float(VecD[n])
        B_diff[n] = float(VecB[n]) - float(VecD[n])
        C_diff[n] = float(VecC[n]) - float(VecD[n])
        
    #iDir = (B_diff[0]*((B_diff[1]*C_diff[2]) - (B_diff[2]*C_diff[1]))) 
    #jDir = -(B_diff[1]*((B_diff[0]*C_diff[2]) - (B_diff[2]*C_diff[0]))) 
    #kDir = (B_diff[2]*((B_diff[0]*C_diff[1]) - (B_diff[1]*C_diff[0]))) 
    
    iDir = B_diff[1]*C_diff[2] - B_diff[2]*C_diff[1]
    jDir = B_diff[2]*C_diff[0] - B_diff[0]*C_diff[2]
    kDir = B_diff[0]*C_diff[1] - B_diff[1]*C_diff[0]

    productAB = [iDir, jDir, kDir]
    
    det = abs(A_diff[0]*productAB[0] + A_diff[1]*productAB[1] + A_diff[2]*productAB[2])
    return det
    
    #Volume_matrix = {}
    #Volume_matrix[0] = []
    #Volume_matrix[1] = []
    #Volume_matrix[2] = []
    #Volume_matrix[3] = []
    #n = 0
    #for index in Tetrahedron_array[j]:
    #    vertex = Tetrahedron_array[j][vertex]
    #    Volume_matrix[0][n] = vertex[0]
    #    Volume_matrix[1][n] = vertex[1]
    #    Volume_matrix[2][n] = vertex[2]
    #    Volume_matrix[3][n] = vertex[3]
        
        
    
    ##det = (Tetrahedron_matrix[j]).det()
    #return det

def circumcenter(j):
    x = {}
    y = {}
    z = {}
    A_matrix = {}
    for i in range(4):
        x[i] = {}
        y[i] = {}
        z[i] = {}
        A_matrix[i] = {}

    v = 0
    for vertex in Tetrahedrons[j]: #range(len(Tetrahedrons[j])):
        vert = int(vertex)
        x[v][0] = pow(float(P_dict[I][vert][2]), 2) + pow(float(P_dict[I][vert][1]), 2) + pow(float(P_dict[I][vert][0]), 2)
        x[v][1] = float(P_dict[I][vert][1])
        x[v][2] = float(P_dict[I][vert][2])
        x[v][3] = 1
        y[v][0] = pow(float(P_dict[I][vert][2]), 2) + pow(float(P_dict[I][vert][1]), 2) + pow(float(P_dict[I][vert][0]), 2)
        y[v][1] = float(P_dict[I][vert][0])
        y[v][2] = float(P_dict[I][vert][2])
        y[v][3] = 1
        z[v][0] = pow(float(P_dict[I][vert][2]), 2) + pow(float(P_dict[I][vert][1]), 2) + pow(float(P_dict[I][vert][0]), 2)
        z[v][1] = float(P_dict[I][vert][0])
        z[v][2] = float(P_dict[I][vert][1])
        z[v][3] = 1
        A_matrix[v][0] = float(P_dict[I][vert][0])
        A_matrix[v][1] = float(P_dict[I][vert][1])
        A_matrix[v][2] = float(P_dict[I][vert][2])
        A_matrix[v][3] = 1
        v = v + 1
                              
    x_det = (x[0][0]*(x[1][1]*(x[2][2]*x[3][3] - x[2][3]*x[3][2]) - x[2][1]*(x[1][2]*x[3][3] - x[1][3]*x[3][2]) + x[3][1]*(x[1][2]*x[2][3] - x[1][3]*x[2][2]))
            - x[1][0]*(x[0][1]*(x[2][2]*x[3][3] - x[2][3]*x[3][2]) - x[2][1]*(x[0][2]*x[3][3] - x[0][3]*x[3][2]) + x[3][1]*(x[0][2]*x[2][3] - x[0][3]*x[2][2]))
            + x[2][0]*(x[0][1]*(x[1][2]*x[3][3] - x[1][3]*x[3][2]) - x[1][1]*(x[0][2]*x[3][3] - x[0][3]*x[3][2]) + x[3][1]*(x[0][2]*x[1][3] - x[0][3]*x[1][2]))
            - x[3][0]*(x[0][1]*(x[1][2]*x[2][3] - x[1][3]*x[2][2]) - x[1][1]*(x[0][2]*x[2][3] - x[0][3]*x[2][2]) + x[2][1]*(x[0][2]*x[1][3] - x[0][3]*x[1][2]))
            )
        
    y_det = (y[0][0]*(y[1][1]*(y[2][2]*y[3][3] - y[2][3]*y[3][2]) - y[2][1]*(y[1][2]*y[3][3] - y[1][3]*y[3][2]) + y[3][1]*(y[1][2]*y[2][3] - y[1][3]*y[2][2]))
            - y[1][0]*(y[0][1]*(y[2][2]*y[3][3] - y[2][3]*y[3][2]) - y[2][1]*(y[0][2]*y[3][3] - y[0][3]*y[3][2]) + y[3][1]*(x[0][2]*x[2][3] - y[0][3]*y[2][2]))
            + y[2][0]*(y[0][1]*(y[1][2]*y[3][3] - y[1][3]*y[3][2]) - y[1][1]*(y[0][2]*y[3][3] - y[0][3]*y[3][2]) + y[3][1]*(x[0][2]*x[1][3] - y[0][3]*y[1][2]))
            - y[3][0]*(y[0][1]*(y[1][2]*y[2][3] - y[1][3]*y[2][2]) - y[1][1]*(y[0][2]*y[2][3] - y[0][3]*y[2][2]) + y[2][1]*(x[0][2]*x[1][3] - y[0][3]*y[1][2]))
            )
    
    z_det = (z[0][0]*(z[1][1]*(z[2][2]*z[3][3] - z[2][3]*z[3][2]) - z[2][1]*(z[1][2]*z[3][3] - z[1][3]*z[3][2]) + z[3][1]*(z[1][2]*z[2][3] - z[1][3]*z[2][2]))
            - z[1][0]*(z[0][1]*(z[2][2]*z[3][3] - z[2][3]*z[3][2]) - z[2][1]*(z[0][2]*z[3][3] - z[0][3]*z[3][2]) + z[3][1]*(z[0][2]*z[2][3] - z[0][3]*z[2][2]))
            + z[2][0]*(z[0][1]*(z[1][2]*z[3][3] - z[1][3]*z[3][2]) - z[1][1]*(z[0][2]*z[3][3] - z[0][3]*z[3][2]) + z[3][1]*(z[0][2]*z[1][3] - z[0][3]*z[1][2]))
            - z[3][0]*(z[0][1]*(z[1][2]*z[2][3] - z[1][3]*z[2][2]) - z[1][1]*(z[0][2]*z[2][3] - z[0][3]*z[2][2]) + z[2][1]*(z[0][2]*z[1][3] - z[0][3]*z[1][2]))
            )
                              
    A_det = (A_matrix[0][0]*(A_matrix[1][1]*(A_matrix[2][2]*A_matrix[3][3] - A_matrix[2][3]*A_matrix[3][2]) - A_matrix[2][1]*(A_matrix[1][2]*A_matrix[3][3] - A_matrix[1][3]*A_matrix[3][2]) + A_matrix[3][1]*(A_matrix[1][2]*A_matrix[2][3] - A_matrix[1][3]*A_matrix[2][2]))
            - A_matrix[1][0]*(A_matrix[0][1]*(A_matrix[2][2]*A_matrix[3][3] - A_matrix[2][3]*A_matrix[3][2]) - A_matrix[2][1]*(A_matrix[0][2]*A_matrix[3][3] - A_matrix[0][3]*A_matrix[3][2]) + A_matrix[3][1]*(A_matrix[0][2]*A_matrix[2][3] - A_matrix[0][3]*A_matrix[2][2]))
            + A_matrix[2][0]*(A_matrix[0][1]*(A_matrix[1][2]*A_matrix[3][3] - A_matrix[1][3]*A_matrix[3][2]) - A_matrix[1][1]*(A_matrix[0][2]*A_matrix[3][3] - A_matrix[0][3]*A_matrix[3][2]) + A_matrix[3][1]*(A_matrix[0][2]*A_matrix[1][3] - A_matrix[0][3]*A_matrix[1][2]))
            - A_matrix[3][0]*(A_matrix[0][1]*(A_matrix[1][2]*A_matrix[2][3] - A_matrix[1][3]*A_matrix[2][2]) - A_matrix[1][1]*(A_matrix[0][2]*A_matrix[2][3] - A_matrix[0][3]*A_matrix[2][2]) + A_matrix[2][1]*(A_matrix[0][2]*A_matrix[1][3] - A_matrix[0][3]*A_matrix[1][2]))
            )

    c = {}
    c[0] = x_det/(2*A_det)
    c[1] = - y_det/(2*A_det)
    c[2] = z_det/(2*A_det)
    return c

A_values = {}

def mapping_A(A):
    A.assemble()
    #N = A.size[0]
    ##A0 = pc.Mat()
    ##A0.create(pc.COMM_WORLD)
    ##A0.setSizes([3*N, 3*N])
    ##A0.setUp()
    ##A0.zeroEntries()
    #for m in range(N):
    #    for n in range(N):
    #        A0[m, n] = A[m, n]
    #for m in range(N, 2*N):
    #    for n in range(N, 2*N):
    #        A0[m, n] = A[m, n]
    #for m in range(2*N, 3*N):
    #    for n in range(2*N, 3*N):
    #        A0[m, n] = A[m, n]
    ##for m in range(N):
    ##    for n in range(N):
    ##        A0[m, n] = A[m, n]
    ##        A0[m, N + n] = A[m, n]
    ##        A0[m, 2*N + n] = A[m, n]
    A0 = pc.Mat()
    A0.create(pc.COMM_WORLD)
    A0.setSizes([3*N, 3*N])
    A0.setUp()
    A0.zeroEntries()
    N_1 = []
    #N_2 = []
    #N_3 = []
    #for n in range(N):
    #    N_1.append(n)
    #    N_2.append(N+n)
    #    N_3.append(2*N+n)
    #print('1')
    #A0.setValuesBlocked(N_1, N_1, A.getValues(N_1, N_1))
    #print('2')
    #A0.setValuesBlocked(N_2, N_2, A.getValues(N_1, N_1))
    #print('3')
    #A0.setValuesBlocked(N_3, N_3, A.getValues(N_1, N_1))
    #print('4')
    #A0.assemble()
    #values = A.getValues(N_1, N_1)
    for i in range(N):#len(values)):
        #for j in range(len(values[i])):
        #for J in range(len(A_values[i])):
        for J in A_values[i]:
        	#j = int(j_str)
        	j = int(J)
        	if j < N:# and j != 0:
        	    if A[i, j] != 0:#values[i][j] != 0:
        		    A0[i, j] = A[i, j]#values[i][j]
        		    A0[N+i, N+j] = A[i, j]#values[i][j]
        		    A0[2*N+i, 2*N+j] = A[i, j]#values[i][j]
    A0.assemble()
    return A0

def mapping_E(g_E):
    #N = g_E.size
    g_E0 = pc.Vec()
    g_E0.create(pc.COMM_WORLD)
    g_E0.setSizes(3*N)
    g_E0.setUp()
    g_E0.zeroEntries()
    for x in range(N):
        g_E0[x] = g_E[x][0]
    for y in range(N):
        g_E0[y+N] = g_E[y][1]
    for z in range(N):
        g_E0[z+2*N] = g_E[z][2]
    return g_E0

def inverse_mapping(h0):
    N = h0.size/3
    h = pc.Mat()
    h.create(pc.COMM_WORLD)
    h.setSizes([N, 3])
    h.setUp()
    for x in range(int(N)):
        h[x, 0] = h0[x]
    for y in range(int(N)):
        h[y, 1] = h0[N+y]
    for z in range(int(N)):
        h[z, 2] = h0[2*N+z]
    h.assemble()
    return h

d = 3 ###number of dimensions
#g_E = pc.Mat()
#g_E.create(pc.COMM_WORLD)
#g_E.setSizes(N, 3)
#g_E.setUp()
time2 = time.time()
g_E = {}
for n in range(N):
    w = star(str(n)) ###have to look at xml file format to create a star function
    Sum = pc.Vec()
    Sum.create(pc.COMM_WORLD)
    Sum.setSizes(3)
    Sum.setUp()
    Sum.zeroEntries()
    for j in range(len(w)):
        c = circumcenter(j) ###have to create a circumcenter function
        T = volume(j)
        c_vec = pc.Vec()
        c_vec.create(pc.COMM_WORLD)
        c_vec.setSizes(3)
        c_vec.setUp()
        for index in range(len(c)):
            c_vec[index] = c[index]
        c_vec.assemble()
        P_vec = pc.Vec()
        P_vec.create(pc.COMM_WORLD)
        P_vec.setSizes(3)
        P_vec.setUp()
        for index in range(len(P_dict[I][n])):
            P_vec[index] = P_dict[I][n][index]
        P_vec.assemble()
        #Sum = Sum + (P[i][n] - c)*T*rho(T) ###what is rho(T)?
        #Sum = Sum + (P[i][n] - c)*T #when rho = 1
        Sum = Sum + (P_vec - c_vec)*T
####g_E[i] = (2/(d+2))*Sum
    g_E[n] = (2/(d+1))*Sum ###when rho = 1

#print('g_E[2345]: ', g_E[2345][0], g_E[2345][1], g_E[2345][2])

I = 0
A_bar.assemble()
A_bar.zeroEntries()
for b in range(2):
#while(g_E > Tol):
    A_bar_dict = {}
    I = I + 1
    P_dict[I] = P_dict[I - 1] ## may need to revise if uses too much memory (rewrite old P's)
    #for m in range(N+Nb):
    #    #if m % 100 == 0:
    #        #print(m)
    #    print(m)
    #    for n in range(N+Nb):
    #        if m != n:
    #            w_m = star(str(m))
    #            w_n = star(str(n))
    #            Sum = 0
    #            for tet in w_m:
    #                if tet in w_n:
    #                    ####T = weighted_det(Tetrahedrons[tet]) ##there is a formula for this in the paper
    #                    T = volume(tet) ###if rho = 1
    #                    Sum = Sum + T
    #            if Sum == 0:
    #                Sum = 1
    #            A_bar[m, n] = (-2/(d*(d+1)))*Sum
    #            #print('here')
    for m in range(N+Nb):
        A_bar_dict[m] = {}
        if m % 100 == 0:
            print(m)
        A_values[m] = {}
        w_m = star(str(m))
        for tet_ in w_m:
            for n in Tetrahedrons[tet_]:
                #A_bar.assemble()
                try:
                    Sum = A_bar_dict[m][n]#A_bar[m, n]
                except:
                    Sum = 0
                for tet in star(n):
                    if tet in w_m:
                        T = volume(tet)
                        Sum = Sum + (-2/(d*(d+1)))*T
                A_bar[m, n] = Sum #(-2/(d*(d+1)))*Sum
                A_bar_dict[m][n] = Sum
                #try:
                # 	A_values[m].append(n)
                #except:
                # 	A_values[m] = []
                # 	A_values[m].append(n)
                A_values[m][n] = 1
    #A_bar.assemble()
    for m in range(N):#+Nb):
        if m % 100 == 0:
            print(m)
        #elif m > 24200:
        #    print(m)
        Sum = 0
        #A_bar.assemble()
        #for n in range(N+Nb):
        for n in A_values[m]:
            if m != int(n):
                Sum = Sum - A_bar_dict[m][n]#A_bar[m, n]
        A_bar[m, m] = Sum
        A_values[m][m] = 1
    print('Done with stage one of ', b)
    A_bar.assemble()
    #for m in range(N):
    #    print(m)
    #    for n in range(N):
    #        A[m,n] = A_bar[m,n]
    #submat_iterator = []
    #for n in range(N):
    #    submat_iterator.append(n)
    #isrc = pc.IS()
    #isrc.createGeneral(submat_iterator)
    #isrc.setPermutation()
    #A = A_bar.createSubMatrix(isrc, isrc)
    A = pc.Mat()
    A.create(pc.COMM_WORLD)
    A.setSizes([N, N])
    A.setUp()
    A.zeroEntries()
    for m in range(N):
        for n in A_values[m]:
            A[m,n] = A_bar_dict[m][str(n)]
    A.assemble()
    print('here')
    A0 = mapping_A(A)
    print('here')
    A0.assemble()
    #if b == 0:
	#    with open('/Users/mvhsan/Papers/mesh_dendritic_spine/A_bar_file_better.txt', 'w+') as file:
	#    	#getvaluesA_bar = []
	#    	#for n in range(N+Nb):
	#    		#getvaluesA_bar.append(n)
	#    	#valuesA_bar = A_bar.getValues(getvaluesA_bar, getvaluesA_bar)
	#    	for i in range(N+Nb):
	#    		for j in range(N+Nb):
	#    			if float(A_bar[i,j]) != 0.0:
	#    				file.write(str(i) + ',' + str(j) + ': ' + str(A_bar[i,j]) + '\n')
	#    with open('/Users/mvhsan/Papers/mesh_dendritic_spine/A0_file_better.txt', 'w+') as file:
	    	#getvaluesA0 = []
	    	#for n in range(3*N):
	    	#	getvaluesA0.append(n)
	    	#valuesA0 = A0.getValues(getvaluesA0, getvaluesA0)
	#    	for i in range(3*N):
	#    		for j in range(3*N):
	#    			if float(A0[i,j]) != 0.0:
	#    				file.write(str(i) + ',' + str(j) + ': ' + str(A0[i,j]) + '\n')
    h0 = pc.Vec()
    h0.create(pc.COMM_WORLD)
    h0.setSizes(3*N)
    h0.setUp()
    print('here')
    g_E0 = mapping_E(g_E)
    print('here')
    PC = pc.PC()
    PC.create(pc.COMM_WORLD)
    PC.setOperators(A=A0)
    PC.setType(pc.PC.Type.HYPRE)
    PC.setUp()
    CG = pc.KSP()
    CG.create(pc.COMM_WORLD)
    CG.setType(pc.KSP.Type.CG)
    CG.setPC(PC)
    print('here')
    CG.solve(g_E0, h0)
    h = inverse_mapping(h0)
    ### Need to compute alpha
    h.assemble()
    if b == 0:
        alpha = 0.00001
    else:
        alpha = .0000000002#35#10#.000125#1e-6
    for e in range(h.size[0]):
        for index in range(h.size[1]):
            P_dict[I][e][index] = float(P_dict[I][e][index]) - float(alpha*h[e, index])

#g_E = pc.Mat()
#g_E.create(pc.COMM_WORLD)
#for n in range(N):
#		w = star(n) ###have to look at xml file format to create a star function
#		Sum = vector(3)
#		Sum.zero_entries()
#		for j in range(len(w)):
#			c = circumcenter(j) ###have to create a circumcenter function
#			T = determinant(j)
#			####Sum = Sum + (P[i][n] - c)*T*rho(T) ###what is rho(T)?
#			Sum = Sum + (P[i][n] - c)*T ##when rho = 1
####g_E[i] = (2/(d+2))*Sum
#		g_E[i] = (2/(d+1))*Sum ###when rho = 1
    if b == 1:
        time2 = time.time() - time2
        print("###########################TIME############\n" + str(time2))
    if b != 1:
        g_E = {}
        for n in range(N):
            w = star(str(n)) ###have to look at xml file format to create a star function
            Sum = pc.Vec()
            Sum.create(pc.COMM_WORLD)
            Sum.setSizes(3)
            Sum.setUp()
            Sum.zeroEntries()
            for j in range(len(w)):
                c = circumcenter(j) ###have to create a circumcenter function
                T = volume(j)
                c_vec = pc.Vec()
                c_vec.create(pc.COMM_WORLD)
                c_vec.setSizes(3)
                c_vec.setUp()
                for index in range(len(c)):
                    c_vec[index] = c[index]
                c_vec.assemble()
                P_vec = pc.Vec()
                P_vec.create(pc.COMM_WORLD)
                P_vec.setSizes(3)
                P_vec.setUp()
                for index in range(len(P_dict[I][n])):
                    P_vec[index] = P_dict[I][n][index]
                P_vec.assemble()
                #Sum = Sum + (P[i][n] - c)*T*rho(T) ###what is rho(T)?
                #Sum = Sum + (P[i][n] - c)*T #when rho = 1
                Sum = Sum + (P_vec - c_vec)*T
    ####g_E[i] = (2/(d+2))*Sum
            g_E[n] = (2/(d+1))*Sum ###when rho = 1
    print('Stage ', b, ' done')

#n = 0
#os.system('cp /Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml /Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_smoothed.xml')

#oldfile = open('tetmeshFromBlender_uninverted.xml')
oldfile = open('tt-sr-mit.tet_mesh_uninverted.xml')
oldfile.seek(0)
#with open('tetmeshFromBlender_smoothed_ODT.xml', 'w+') as newfile:
with open('tt-sr-mit.tet_mesh_smoothed_ODT.xml', 'w+') as newfile:
    for line in oldfile:
        if '<vertex' in line:
            index = (re.findall('index="([0-9]*)', line)[0])
            newline = '      <vertex index="{0}" x="{1}" y="{2}" z="{3}" />\n'.format(index, P_dict[I][int(Intermediary[index])][0], P_dict[I][int(Intermediary[index])][1], P_dict[I][int(Intermediary[index])][2])
            #newline = '<vertex index="{0}" x="{1}" y="{2}" z="{3}" />'.format(index, P_dict[I][int(Intermediary[int(index)])][0], P_dict[I][int(Intermediary[int(index)])][1], P_dict[I][int(Intermediary[int(index)])][2])		### I should at this point be the last iteration value
            #n ++
        else:
            newline = line
        newfile.write(newline)
oldfile.close()

#mesh = Mesh('tetmeshFromBlender_smoothed_ODT.xml')
#vtkfile = File('tetmeshFromBlender_smoothed_ODT.pvd')
mesh = Mesh('tt-sr-mit.tet_mesh_smoothed_ODT.xml')
vtkfile = File('tt-sr-mit.tet_mesh_smoothed_ODT.pvd')
vtkfile << mesh
