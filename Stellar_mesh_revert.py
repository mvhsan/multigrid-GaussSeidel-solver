import re
from dolfin import *

Vertices = {}
with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted.node') as file:
	n = 0
	for line in file:
		if n == 0:
			v_size = re.findall('([0-9]*) ', line)[0]
			n = n + 1
		else:
			index = re.findall('([0-9]*) ', line)[0]
			x = re.findall(' ((?:\-|)(?:[0-9]+\.[0-9]+|[0-9]+))', line)[0]
			y = re.findall(' ((?:\-|)(?:[0-9]+\.[0-9]+|[0-9]+))', line)[1]
			z = re.findall(' ((?:\-|)(?:[0-9]+\.[0-9]+|[0-9]+))', line)[2]
			Vertices[index] = [x, y, z]
print('Projected vertices size: ', v_size)
print('Actual vertices size: ', len(Vertices))


Tetrahedrons = {}
with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted.ele') as file:
	n = 0
	for line in file:
		if n == 0:
			t_size = re.findall('([0-9]*) ', line)[0]
			n = n + 1
		else:
			index = re.findall('([0-9]*) ', line)[0]
			v0 = re.findall(' ([0-9]+)', line)[0]
			v1 = re.findall(' ([0-9]+)', line)[1]
			v2 = re.findall(' ([0-9]+)', line)[2]
			v3 = re.findall(' ([0-9]+)', line)[3]
			#Tetrahedrons[index] = [v1, v2, v3, v0]
			Tetrahedrons[index] = [v0, v1, v2, v3]

print('Projected tetrahedra size: ', t_size)
print('Actual tetrahedra size: ', len(Tetrahedrons))

oldfile = open('/Users/mvhsan/GAMer/2018_DSM_GAMer/checkpoints/tetmeshFromBlender.xml')
with open('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_testing.xml', 'w+') as file:
	#file.write('<?xml version="1.0"?>\n<dolfin xmlns:dolfin="http://fenicsproject.org">\n  <mesh celltype="tetrahedron" dim="3">\n')
	#file.write('    <vertices size="' + v_size + '">\n')
	#for vertex in Vertices:
#		line = ('      <vertex index="' + vertex + '" x="' + Vertices[vertex][0] + '" y="' + Vertices[vertex][1] + '" z="' + Vertices[vertex][2] + '" />\n')
#		file.write(line)
#	file.write('    </vertices>\n')
#	file.write('    <cells size="' + t_size + '">\n')
#	for tetrahedron in Tetrahedrons:
#		line = ('      <tetrahedron index="' + tetrahedron + '" v0="' + Tetrahedrons[tetrahedron][0] + '" v1="' + Tetrahedrons[tetrahedron][1] + '" v2="' + Tetrahedrons[tetrahedron][2] + '" v3="' + Tetrahedrons[tetrahedron][3] + '" />\n')
#		file.write(line)
#	file.write('    </cells>\n')
#	file.write('    <data />\n')
#	file.write('  </mesh>\n')
#	file.write('</dolfin>')
	tetra = 0
	vert = 0
	for line in oldfile:
		if '<vertex ' in line:
			newline = ('      <vertex index="' + str(vert) + '" x="' + Vertices[str(vert)][0] + '" y="' + Vertices[str(vert)][1] + '" z="' + Vertices[str(vert)][2] + '" />\n')
			vert = vert + 1
		elif '<tetrahedron ' in line:
			newline = ('      <tetrahedron index="' + str(tetra) + '" v0="' + Tetrahedrons[str(tetra)][0] + '" v1="' + Tetrahedrons[str(tetra)][1] + '" v2="' + Tetrahedrons[str(tetra)][2] + '" v3="' + Tetrahedrons[str(tetra)][3] + '" />\n')
			tetra = tetra + 1
		else:
			newline = line
		file.write(newline)
oldfile.close()

mesh = Mesh('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_testing.xml')

vtkfile = File('/Users/mvhsan/Papers/mesh_dendritic_spine/tetmeshFromBlender_uninverted_testing.pvd')
vtkfile << mesh
