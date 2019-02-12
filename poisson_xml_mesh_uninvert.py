######################Convert .xml file to .node and .ele files################

import re

Tetrahedrons_ = {}
#with open('/Users/mvhsan/Papers/dendritic_spine_meshing/Stellar_files/testmesh.xml') as file:
#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.xml') as file:
#with open('tetmeshFromBlender.xml') as file:
with open('circlemesh.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.xml') as file:
    for line in file:
        if '<triangle ' in line:
            a_ = re.findall('index="([0-9]*)', line)
            a = a_[0]  ##or whatever the syntax is
            v0_ = re.findall('v0="([0-9]*)', line)  ##need to fix syntax here
            v1_ = re.findall('v1="([0-9]*)', line)
            v2_ = re.findall('v2="([0-9]*)', line)
            #v3_ = re.findall('v3="([0-9]*)', line)
            v0 = v0_[0]
            v1 = v1_[0]
            v2 = v2_[0]
            #v3 = v3_[0]
            #Tetrahedrons_[a] = [v3, v0, v1, v2]
            Tetrahedrons_[a] = [v0, v1, v2]#, #v3]

Original_vertices = {}
#n = 0
N = 0
#with open('/Users/mvhsan/Papers/dendritic_spine_meshing/Stellar_files/testmesh.xml') as file:
#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.xml') as file:
#with open('tetmeshFromBlender.xml') as file:
with open('circlemesh.xml') as file:
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
            #z = (re.findall(' z="((?:\-|)(?:[0-9]*\.[0-9]*|[0-9]*))', line))[0]
            try:
                xE = (re.findall('e((?:\+|\-)[0-9]*)', line))[0]
                yE = (re.findall('e((?:\+|\-)[0-9]*)', line))[1]
                #zE = (re.findall('e((?:\+|\-)[0-9]*)', line))[2]
            except:
                xE = 0
                yE = 0
                #zE = 0
            #print(x + xE)
            x_num = float(x)*pow(10, float(xE))
            y_num = float(y)*pow(10, float(yE))
            #z_num = float(z)*pow(10, float(zE))
            #print(x_num)
            #print(str(x_num))
            #print(x_num)
            #print(re.findall(' x="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))
            #x = (re.findall(' x="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))[0]
            #y = (re.findall(' y="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))[0]
            #z = (re.findall(' z="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))[0]
            Original_vertices[index] = [str(x_num), str(y_num)]#, str(z_num)]
            #print(Original_vertices[index])

#with open('/Users/mvhsan/Papers/dendritic_spine_meshing/Stellar_files/testmesh.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.xml') as file:
#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.xml') as file:
#with open('tetmeshFromBlender.xml') as file:
with open('circlemesh.xml') as file:
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
    #D = 0
    Boundary_bool = {}
    for line in file:
        if itr > 0:
            #D = D + 1
            if '<value' in line:
                #D = D + 1
                if 'value="1"' not in line:
                #if 'value="21"' in line:
                #if 'value="1"' in line:
                    entity = (re.findall('local_entity="([0-9]*)', line))[0]
                    tet = (re.findall('cell_index="([0-9]*)', line))[0]
                    for vert in Tetrahedrons_[tet]:
                        if vert != entity:
                            if vert not in Boundary_bool:
                                Boundary_bool[str(vert)] = 1
            if 'mesh_value_collection' and 'dim="3"' in line:
                #print('el')
                break

for vert in range(N):
    if str(vert) not in Boundary_bool:
        Boundary_bool[str(vert)] = 0

# def volume(j):
#     #p = pc.PC()
#     #p.create(pc.COMM_WORLD)
#     #p.setType(pc.PC.Type.LU)
#     #p.setOperators(A=Tetrahedron_matrix[j])
#     #p.setUp()

#     #k = pc.KSP()
#     #k.create(pc.COMM_WORLD)
#     #k.setPC(p)

#     VecA = Tetrahedron_array[j][0]
#     VecB = Tetrahedron_array[j][1]
#     VecC = Tetrahedron_array[j][2]
#     VecD = Tetrahedron_array[j][3]

#     A_diff = {}
#     B_diff = {}
#     C_diff = {}
#     for n in range(3):
#         A_diff[n] = float(VecA[n]) - float(VecD[n])
#         B_diff[n] = float(VecB[n]) - float(VecD[n])
#         C_diff[n] = float(VecC[n]) - float(VecD[n])

#     #iDir = (B_diff[0]*((B_diff[1]*C_diff[2]) - (B_diff[2]*C_diff[1])))
#     #jDir = -(B_diff[1]*((B_diff[0]*C_diff[2]) - (B_diff[2]*C_diff[0])))
#     #kDir = (B_diff[2]*((B_diff[0]*C_diff[1]) - (B_diff[1]*C_diff[0])))

#     iDir = B_diff[1]*C_diff[2] - B_diff[2]*C_diff[1]
#     jDir = B_diff[2]*C_diff[0] - B_diff[0]*C_diff[2]
#     kDir = B_diff[0]*C_diff[1] - B_diff[1]*C_diff[0]
#     productAB = [iDir, jDir, kDir]

#     det = (A_diff[0]*productAB[0] + A_diff[1]*productAB[1] + A_diff[2]*productAB[2])/6
#     return det

# Tetrahedron_array = {}
# #t_n = 0
# for tetrahedron in Tetrahedrons_:
#     Tet = {}
#     m = 0
#     for vertex in Tetrahedrons_[tetrahedron]:
#         Tet[m] = Original_vertices[vertex]
#         m = m + 1
#     Tetrahedron_array[tetrahedron] = Tet
#     #t_n = t_n + 1

# Tetrahedrons = {}
# for tetrahedron in Tetrahedrons_:
#     v0 = Tetrahedrons_[tetrahedron][0]
#     v1 = Tetrahedrons_[tetrahedron][1]
#     v2 = Tetrahedrons_[tetrahedron][2]
    #v3 = Tetrahedrons_[tetrahedron][3]
    #if volume(tetrahedron) < 0:
    #    Tetrahedrons[tetrahedron] = [v3, v0, v1, v2]
    #else:
    #    Tetrahedrons[tetrahedron] = [v0, v1, v2, v3]

#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.node', 'w+') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.node', 'w+') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.node', 'w+') as file:
#print(Boundary_bool)
#print(Tetrahedrons_)

#oldfile = open('tetmeshFromBlender.xml')
oldfile = open('circlemesh.xml')
oldfile.seek(0)
#with open('tetmeshFromBlender_uninverted.xml', 'w+') as newfile:
with open('circlemesh_uninverted.xml', 'w+') as newfile:
    for line in oldfile:
        if '<vertex' in line:
            index = (re.findall('index="([0-9]*)', line)[0])
            newline = '      <vertex index="{0}" x="{1}" y="{2}" />\n'.format(index, Original_vertices[index][0], Original_vertices[index][1])
        elif '<triangle' in line:
            index = (re.findall('index="([0-9]*)', line)[0])
            newline = '      <triangle index="{0}" v0="{1}" v1="{2}" v2="{3}" />\n'.format(index, Tetrahedrons_[index][0], Tetrahedrons_[index][1], Tetrahedrons_[index][2])
            #newline = '<vertex index="{0}" x="{1}" y="{2}" z="{3}" />'.format(index, P_dict[I][int(Intermediary[int(index)])][0], P_dict[I][int(Intermediary[int(index)])][1], P_dict[I][int(Intermediary[int(index)])][2])        ### I should at this point be the last iteration value
            #n ++
        else:
            newline = line
        newfile.write(newline)
oldfile.close()
