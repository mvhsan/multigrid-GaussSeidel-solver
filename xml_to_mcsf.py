#convert to mc-compatible file


import re

Tetrahedrons_ = {}
#with open('/Users/mvhsan/Papers/dendritic_spine_meshing/Stellar_files/testmesh.xml') as file:
#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
#with open('tt-sr-mit.tet_mesh_smoothed_ODT.xml') as file:
with open('circlemesh_uninverted.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.xml') as file:
    for line in file:
        if '<tetrahedron ' in line:
            a_ = re.findall('index="([0-9]*)', line)
            a = a_[0]  ##or whatever the syntax is
            v0_ = re.findall('v0="([0-9]*)', line)  ##need to fix syntax here
            v1_ = re.findall('v1="([0-9]*)', line)
            v2_ = re.findall('v2="([0-9]*)', line)
            v3_ = re.findall('v3="([0-9]*)', line)
            v0 = v0_[0]
            v1 = v1_[0]
            v2 = v2_[0]
            v3 = v3_[0]
            #Tetrahedrons_[a] = [v3, v0, v1, v2]
            Tetrahedrons_[a] = [v0, v1, v2, v3]

Original_vertices = {}
#n = 0
N = 0
#with open('/Users/mvhsan/Papers/dendritic_spine_meshing/Stellar_files/testmesh.xml') as file:
#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
with open('circlemesh_uninverted.xml') as file:
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
            try:
                xE = (re.findall('e((?:\+|\-)[0-9]*)', line))[0]
                yE = (re.findall('e((?:\+|\-)[0-9]*)', line))[1]
                zE = (re.findall('e((?:\+|\-)[0-9]*)', line))[2]
            except:
                xE = 0
                yE = 0
                zE = 0
            #print(x + xE)
            x_num = float(x)*pow(10, float(xE))
            y_num = float(y)*pow(10, float(yE))
            z_num = float(z)*pow(10, float(zE))
            #print(x_num)
            #print(str(x_num))
            #print(x_num)
            #print(re.findall(' x="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))
            #x = (re.findall(' x="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))[0]
            #y = (re.findall(' y="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))[0]
            #z = (re.findall(' z="((?:\-|)(?:[0-9]*\.[0-9]*e(?:\-|\+)[0-9]*|[0-9]*))', line))[0]
            Original_vertices[index] = [str(x_num), str(y_num), str(z_num)]
            #print(Original_vertices[index])

#with open('/Users/mvhsan/Papers/dendritic_spine_meshing/Stellar_files/testmesh.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/testing.xml') as file:
#with open('/Users/mvhsan/tetmeshFromBlender_smoothed.xml') as file:
#with open('/Users/mvhsan/Stellar_1.0/cubemesh_9-46656.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
with open('circlemesh_uninverted.xml') as file:
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

#with open('tetmeshFromBlender_smoothed_ODT.m', 'w+') as file:
with open('circlemesh_uninverted.m', 'w+') as file:
	file.write('mcsf_begin=1;\n\n')
	file.write('      dim=3;\n      dimii=3;\n      vertices=' + str(N) + ';\n      simplices=' + str(len(Tetrahedrons_)) + ';\n\n')
	file.write('vert=[\n')
	for vert in Original_vertices:
		file.write(vert + '    0    ' + Original_vertices[vert][0] + '  ' + Original_vertices[vert][1] + '  ' + Original_vertices[vert][2] + '\n')
	file.write('];\n\n')
	file.write('simp=[\n')
	for tetra in Tetrahedrons_:
		file.write(tetra + '      0   0   ' + str(Boundary_bool[Tetrahedrons_[tetra][0]]) + '  ' + str(Boundary_bool[Tetrahedrons_[tetra][1]]) + '  ' + str(Boundary_bool[Tetrahedrons_[tetra][2]]) + '  ' + str(Boundary_bool[Tetrahedrons_[tetra][3]]) + '  ' + Tetrahedrons_[tetra][0] + '  ' + Tetrahedrons_[tetra][1] + '  ' + Tetrahedrons_[tetra][2] + '  ' + Tetrahedrons_[tetra][3] + '\n')
	file.write('];\n\n')
	file.write('mcsf_end=1;')