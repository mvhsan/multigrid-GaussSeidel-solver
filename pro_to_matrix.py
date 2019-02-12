#make petsc matrix from prolongation .m file

import petsc4py.PETSc as pc
import re

with open('pro.m') as file:
	for line in file:
		if 'dim=' in line:
			dim1 = re.findall('([0-9]+)x[0-9]+', line)[0]
			dim2 = re.findall('[0-9]+x([0-9]+)', line)[0]
			break
	A = pc.Mat()
	A.create(pc.COMM_WORLD)
	A.setSizes([int(dim1), int(dim2)])
	A.setUp()
	for line in file:
		if 'e+' or 'e-' in line:
			i = re.findall('  ([0-9]+)  ', line)[0]
			j = re.findall('  ([0-9]+)  ', line)[1]
			val = re.findall('((?:\-|)[0-9]+\.[0-9]+)', line)[0]
			power = re.findall('e(?:\-|)([0-9]+)', line)[0]
			A[i, j] = val * pow(10, power)

A.assemble()
