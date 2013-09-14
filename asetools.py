#!/usr/bin/python

"""Extra functions to compliment the ASE package. All functions work on an ASE Atoms object. """
import numpy
import matplotlib.pyplot as pyplot

def xyz(coordinates,atoms):
	"""Converts an array of fractional coordinates to xyz coordinates for a given cell """
	coord_list =[]
	for coordinate in coordinates:
		fcoord = numpy.dot(coordinate,atoms.get_cell())
		coord_list.append(fcoord)
	xyz_coordinates = numpy.vstack(coord_list)
	return xyz_coordinates

def fractional(coordinates,atoms):
	"""Converts an array of xyz coordinates to fractional coordinates for a given Atoms object """
	#import pdb;pdb.set_trace()
	coord_list =[]
	for coordinate in coordinates:
		fcoord = numpy.dot(coordinate,atoms.get_reciprocal_cell().T)
		coord_list.append(fcoord)
	fractional_coordinates = numpy.vstack(coord_list)
	return fractional_coordinates

def magnitude(array):
	""" This returns the magnitude of vector properties from ASE. It calculates the magnitude or vector norm of 1D slices in the given array. """
	mag = numpy.apply_along_axis(numpy.linalg.norm,1,array)
	return mag

def cellparam(cell):
	""" Returns a 6-element vector that represents the cell parameters, (a,b,c,alpha,beta,gamma).
	Required input is an Atoms object from ASE."""
	abc  = numpy.sqrt(numpy.sum(cell.get_cell()**2,axis=1))
	alpha = numpy.float64( numpy.degrees ( numpy.arccos( numpy.dot( cell.get_cell()[2:3], cell.get_cell()[1:2].T) / (( abc[1]*abc[2]))) ))
	beta = numpy.float64( numpy.degrees ( numpy.arccos( numpy.dot( cell.get_cell()[0:1], cell.get_cell()[2:3].T) / (( abc[0]*abc[2]))) ))
	gamma = numpy.float64( numpy.degrees ( numpy.arccos( numpy.dot( cell.get_cell()[0:1], cell.get_cell()[1:2].T) / (( abc[0]*abc[1]))) ))
	return [abc[0], abc[1], abc[2], alpha, beta, gamma]

def clustercentre(defects):
	"""Finds the centre of a defect cluster """
	raise NotImplementedError("Possible future function but yet to be written.")


def scatter_plot(title,x_title,x,x_units,y_title,y,y_units):
	""" plots two arrays against each other """
	# Check arrays are of the same length.
	if (x.shape != y.shape):
		raise Exception('x and y arrays are of different shapes')
	fignum = random.randint(0,sys.maxint)
	pyplot.figure(fignum)
	### Plotting arrangements ###
	pyplot.scatter(x,y)
	pyplot.title(title)
	pyplot.xlabel(xtitle)
	pyplot.ylabel(ytitle)
	raise NotImplementedError("Only partially written.")
	return fignum

def get_chg_array(file):
	""" Returns the charge densities on the fine FFT-grid in a numpy array."""
	from ase.calculators.vasp import VaspChargeDensity
	pot_array = numpy.vstack((VaspChargeDensity(file)).chg)
	return pot_array

def get_locpot_array(file):
	""" Returns the charge densities on the fine FFT-grid in a numpy array."""
	from ase.calculators.vasp import VaspChargeDensity
	from ase.calculators.vasp import Vasp
 	calc = Vasp(restart=True)
 	cell = calc.get_atoms()
 	vol = cell.get_volume()
	pot_array = numpy.vstack((VaspChargeDensity(file)).chg)
	numpy.divide(pot_array,vol)
	return pot_array

def get_average_local_potential(file):
	""" Returns the average potential on the fine FFT-grid as given in the file specified."""
	pot_array = get_locpot_array(file)
	average = numpy.mean(pot_array)
	return average

def get_distances_to_defect(number):
	""" Returns the distances to the atom number "number" """
	calc = ase.calculators.vasp.Vasp(restart=1)
	atoms = calc.get_atoms()
	aseroutines.output_distances_to_atom(atoms,number)

def get_distances_to_point(a,b,c):
	""" Returns the distances between each atom and the coordinate "a,b,c" """
	calc = ase.calculators.vasp.Vasp(restart=1)
	atoms = calc.get_atoms()
	aseroutines.output_distances_to_point(atoms,a,b,c)

def show(matrix):
	# Print out matrix
	for col in matrix:
		print col

def zero(m,n):
	# Create zero matrix
	new_matrix = [[0 for row in range(n)] for col in range(m)]
	return new_matrix

def mult(matrix1,matrix2):
	# Matrix multiplication
	if len(matrix1[0]) != len(matrix2):
		# Check matrix dimensions
		print 'Matrices must be m*n and n*p to multiply!'
	else:
		# Multiply if correct dimensions
		new_matrix = zero(len(matrix1),len(matrix2[0]))
		for i in range(len(matrix1)):
			for j in range(len(matrix2[0])):
				for k in range(len(matrix2)):
					new_matrix[i][j] += matrix1[i][k]*matrix2[k][j]
		return new_matrix

def symmetrize(matrix):
    return matrix + matrix.T - numpy.diag(matrix.diagonal())


def calc_defect_strain_matrix(Latt):
	from ase.calculators.vasp import Vasp
	# defect lattice
	calc = Vasp(restart=1)
	atoms = calc.get_atoms()
	LattDef = numpy.array(atoms.get_cell())
	dstrain = conc*numpy.dot(numpy.subtract(LattDef,Latt.T),numpy.linalg.inv(Latt.T))
	dstrain=asetools.symmetrize(dstrain)
	return dstrain

def tensorise(matrix):
	tensor = numpy.array([matrix[0][0],matrix[1][1],matrix[2][2],2*matrix[1][2],2*matrix[0][2],2*matrix[0][1]])
	return tensor


def xmgrace_reverse_polarity(filename):
	""" Reverses the polarity of spin-down plots in the xmgrace file supplied """
	import re
	from StringIO import StringIO
	import os
	f = open(filename)
	lines = f.read().splitlines()
	f.close
	os.rename(os.path.realpath(filename), os.path.realpath(filename)+".bak")
	new_lines = []
	for line in lines:
		if re.match("\s*[-0-9.]\s*[-0-9.]", line):
			numbers = line.split()
			numbers[1] = float(numbers[1])*-1.0
			new_line = '{:>22}{:>19}'.format(numbers[0],numbers[1])
			new_lines.append(new_line)
		else:
			new_lines.append(line)
	outfile = filename
	fout = open(outfile,'w')
	for line in new_lines:
		print>>fout, line
	fout.close

