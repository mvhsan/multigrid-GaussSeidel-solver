#############################measure the quality of a mesh################
### if given the Tetrahedron and P objects that I constructed in the ODT code

import re
import math
import numpy
import statistics

inverted = []

Tetrahedrons_ = {}
#with open('/Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
#with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_smoothed.xml') as file:
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
            Tetrahedrons_[a] = [v0, v1, v2, v3]

Original_vertices = {}
#n = 0
N = 0
#with open('/Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
#with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_smoothed.xml') as file:
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
#with open('/Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
#with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_smoothed.xml') as file:
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

Tetrahedrons = {}
#with open('/Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml') as file:
#with open('tetmeshFromBlender_smoothed_ODT.xml') as file:
with open('tt-sr-mit.tet_mesh_uninverted.xml') as file:
#with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_smoothed.xml') as file:
    for line in file:
        if '<tetrahedron' in line:
            a_ = re.findall('index="([0-9]*)', line)
            a = int(a_[0])  ##or whatever the syntax is
            v0_ = re.findall('v0="([0-9]*)', line)  ##need to fix syntax here
            v1_ = re.findall('v1="([0-9]*)', line)
            v2_ = re.findall('v2="([0-9]*)', line)
            v3_ = re.findall('v3="([0-9]*)', line)
            #print(line)
            #print(a_[0] + ' ' + v0_[0] + ' ' + v1_[0] + ' ' + v2_[0] + ' ' + v3_[0])
            v0 = int(v0_[0])
            v1 = int(v1_[0])
            v2 = int(v2_[0])
            v3 = int(v3_[0])
            #Tetrahedrons[a] = [Reverse_intermediary[v0], Reverse_intermediary[v1], Reverse_intermediary[v2], Reverse_intermediary[v3]]
            Tetrahedrons[a] = [Intermediary[str(v0)], Intermediary[str(v1)], Intermediary[str(v2)], Intermediary[str(v3)]]

print(len(Tetrahedrons))

I = 0

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

	det = (A_diff[0]*productAB[0] + A_diff[1]*productAB[1] + A_diff[2]*productAB[2])/6
	if det < 0:
		inverted.append(j)
	return abs(det)

def dist(vert1, vert2):
	x = float(vert1[0]) - float(vert2[0])
	y = float(vert1[1]) - float(vert2[1])
	z = float(vert1[2]) - float(vert2[2])
	
	distance = math.sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
	return distance

def Area(vert1, vert2, vert3):
	AB_mag = dist(vert1, vert2)
	AC_mag = dist(vert1, vert3)
	
	AB = {}
	AC = {}

	AB[0] = float(vert1[0]) - float(vert2[0])
	AB[1] = float(vert1[1]) - float(vert2[1])
	AB[2] = float(vert1[2]) - float(vert2[2])
	AC[0] = float(vert1[0]) - float(vert3[0])
	AC[1] = float(vert1[1]) - float(vert3[1])
	AC[2] = float(vert1[2]) - float(vert3[2])

	ABdotAC = AB[0]*AC[0] + AB[1]*AC[1] + AB[2]*AC[2]
	theta = numpy.arccos(ABdotAC/(AB_mag*AC_mag))

	area = 0.5*AB_mag*AC_mag*math.sin(theta)	
	return area


def quality_measure(j):
	VecA = Tetrahedron_array[j][0]
	VecB = Tetrahedron_array[j][1]
	VecC = Tetrahedron_array[j][2]
	VecD = Tetrahedron_array[j][3]

	A1 = Area(VecA, VecB, VecC)
	A2 = Area(VecA, VecC, VecD)
	A3 = Area(VecA, VecB, VecD)
	A4 = Area(VecB, VecC, VecD)

	V = volume(j)
	
	inradius = (3*V/(A1+A2+A3+A4))

	a_s = dist(VecA, VecB)
	b_s = dist(VecA, VecC)
	c_s = dist(VecA, VecD)
	A_s = dist(VecB, VecC)
	B_s = dist(VecB, VecD)
	C_s = dist(VecC, VecD)

	product = (a_s*A_s + b_s*B_s + c_s*C_s)*(a_s*A_s + b_s*B_s - c_s*C_s)*(a_s*A_s - b_s*B_s + c_s*C_s)*(-a_s*A_s + b_s*B_s + c_s*C_s)
	#circumradius = math.sqrt(product)/(24*V)
	circumradius = math.sqrt(abs(product))/(24*V)
	#if product < 0:
	#	inverted.append(j)

	qual = (inradius/circumradius)
	#qual = (circumradius/inradius)
	return qual

qual_list = []
qual_ave = 0
qual_min = 1000
for j in Tetrahedron_array:
	quality = quality_measure(j)
	qual_ave = qual_ave + quality
	qual_list.append(quality)
	if quality < qual_min:
		qual_min = quality
qual_ave = qual_ave/len(Tetrahedron_array)
qual_med = statistics.median(qual_list)

print('Average: ', qual_ave)
print('Minimum: ', qual_min)
print('Median: ', qual_med)
print('Inverted: ', len(inverted))
