#!/usr/bin/python

"""Routines to manipulate VASP simulations with ASE. """
import numpy
import shutil
import os
import ase
import asetools
from tempfile import mkstemp

def setup_dvol_corrections_evq(self):
 from ase.calculators.vasp import Vasp
 calc = Vasp(restart=True)
 atoms = calc.get_atoms()
 cell = atoms.get_cell()
 #Create small and large lattice matrices
 small = cell*0.99
 large = cell*1.01
 # Setup the large and small POSCAR files
 atoms.set_cell(small, scale_atoms=True)
 ase.io.write('POSCAR99pc',atoms,format='vasp', direct=1)
 atoms.set_cell(large, scale_atoms=True)
 ase.io.write('POSCAR101pc',atoms,format='vasp', direct=1
 )
 # copy the basic VASP files to the new directories
 directories = ['0x99pc','1x00pc','1x01pc']
 files = ['KPOINTS','POTCAR','WAVECAR','CHGCAR']
 for dir in directories:
  if not os.path.exists(dir):
   os.makedirs(dir)
  for file in files:
   if os.path.exists(file):
    shutil.copy(file,dir)
 # copy the new POSCAR files to the new directories
 shutil.copy('POSCAR','1x00pc')
 shutil.copy('POSCAR99pc','0x99pc/POSCAR')
 shutil.copy('POSCAR101pc','1x01pc/POSCAR')
 # Edit INCAR to create a single point calculations, and save in new directories
 EVQ_INCAR = 'INCAR_EVQ'
 bad_words = ['NSW', 'LVTOT', 'NELM ', 'ISTART', 'ICHARG']
 with open('INCAR') as oldfile, open(EVQ_INCAR, 'w') as newfile:
  for line in oldfile:
   if not any(bad_word in line for bad_word in bad_words):
    newfile.write(line)
  newfile.write('ISTART = 0\nICHARGE = 2\nNSW  = 0\nLVTOT  = .TRUE.\nNELM  = 60\n')
 # Save the new INCAR in directories
 shutil.copy(EVQ_INCAR,'1x00pc/INCAR')
 shutil.copy(EVQ_INCAR,'0x99pc/INCAR')
 shutil.copy(EVQ_INCAR,'1x01pc/INCAR')



def edit_LOCPOT_file(filename):
# Removes a spurious set of lines between potential tables for major and minor spins so ASE will read it correctly
	from StringIO import StringIO
	f = open(filename)
	lines = f.read().splitlines()
	n_atoms = numpy.sum(numpy.genfromtxt(StringIO(lines[6]),dtype=(int)))
	n_points = numpy.product(numpy.genfromtxt(StringIO(lines[9+n_atoms]),dtype=(int)))
	f.close
	#get start line of the potential table for majority spin
	start_line = 10+n_atoms

	#get lastline of the potential table for majority spin
	if numpy.remainder(n_points,5) == 0:
		last_line = 9+n_atoms+n_points/5
	elif numpy.remainder(n_points,5) != 0:
		last_line = 9+n_atoms+n_points/5+1
	#print  "%d %s" % (last_line+1, lines[last_line+1])
	if numpy.remainder(n_atoms,5) == 0:
		n_atom_table_lines = n_atoms/5
	elif numpy.remainder(n_atoms,5) != 0:
		n_atom_table_lines = n_atoms/5+1
	last_atom_table_line = n_atom_table_lines+last_line
	#print  "%d %s" % (last_atom_table_line, lines[last_atom_table_line])
	del lines[last_line+1:last_atom_table_line+1]
	outfile = 'editted_LOCPOT_file'
	fout = open(outfile,'w')
	for line in lines:
		print>>fout, line
	fout.close

def get_spin_density_file(filename):
# edits a CHGCAR file to have only the spin difference table, "spin_up - spin_down"
	import re
	from StringIO import StringIO
	f = open(filename)
	lines = f.read().splitlines()
	n_atoms = numpy.sum(numpy.genfromtxt(StringIO(lines[6]),dtype=(int)))
	n_points = numpy.product(numpy.genfromtxt(StringIO(lines[9+n_atoms]),dtype=(int)))
	f.close
	# Remove the first spin table
	#get start line of the potential table for majority spin
	start_line = 10+n_atoms
	#get lastline of the potential table for majority spin
	if numpy.remainder(n_points,5) == 0:
		last_line = 9+n_atoms+n_points/5
	elif numpy.remainder(n_points,5) != 0:
		last_line = 9+n_atoms+n_points/5+1
	del lines[start_line:last_line]

	# delete lines until you next match the "number of grid points line"
	finished = 0;
	count = 0;
	while finished != 1:
		l_match = re.match(lines[9+n_atoms],lines[9+n_atoms+1],0)
		if (l_match):
			finished = 1
		else:
			del lines[9+n_atoms+1]

	del lines[9+n_atoms+1]
	outfile = 'CHGCAR-spin_density'
	fout = open(outfile,'w')
	for line in lines:
		print>>fout, line
	fout.close


def get_EVQ_data():
	data = []
	heading = os.path.realpath('INCAR')
	data.append(heading)
	columns = '{:^18}{:^18}'.format('Volume, A','potential, eV')
	data.append(columns)
	os.chdir('./0x99pc')
	shutil.move('./OUTPUT/LOCPOT','./LOCPOT')
	shutil.move('./OUTPUT/CONTCAR','./CONTCAR')
	shutil.move('./OUTPUT/OUTCAR','./OUTCAR')
	calc = ase.calculators.vasp.Vasp(restart=1)
	vol1 = calc.atoms.get_volume()
	edit_LOCPOT_file('LOCPOT')
	loc = ase.calculators.vasp.VaspChargeDensity('editted_LOCPOT_file')
	avr_pot1 = numpy.mean([loc.chg,loc.chgdiff])*vol1
	small = '{:^18E}{:^18E}'.format(vol1, avr_pot1)
	data.append(small)
	os.chdir('../1x00pc')
	shutil.move('./OUTPUT/LOCPOT','./LOCPOT')
	shutil.move('./OUTPUT/CONTCAR','./CONTCAR')
	shutil.move('./OUTPUT/OUTCAR','./OUTCAR')
	calc = ase.calculators.vasp.Vasp(restart=1)
	vol2 = calc.atoms.get_volume()
	edit_LOCPOT_file('LOCPOT')
	loc = ase.calculators.vasp.VaspChargeDensity('editted_LOCPOT_file')
	avr_pot2 = numpy.mean([loc.chg,loc.chgdiff])*vol2
	standard = '{:^18E}{:^18E}'.format(vol2, avr_pot2)
	data.append(standard)
	os.chdir('../1x01pc')
	shutil.move('./OUTPUT/LOCPOT','./LOCPOT')
	shutil.move('./OUTPUT/CONTCAR','./CONTCAR')
	shutil.move('./OUTPUT/OUTCAR','./OUTCAR')
	calc = ase.calculators.vasp.Vasp(restart=1)
	vol3 = calc.atoms.get_volume()
	edit_LOCPOT_file('LOCPOT')
	loc = ase.calculators.vasp.VaspChargeDensity('editted_LOCPOT_file')
	avr_pot3 = numpy.mean([loc.chg,loc.chgdiff])*vol3
	large = '{:^18E}{:^18E}'.format(vol3, avr_pot3)
	data.append(large)
	os.chdir('../')
	x = numpy.array([vol1,vol2,vol3])
	y = numpy.array([avr_pot1,avr_pot2,avr_pot3])
	fit = numpy.polyfit(x, y, 1)
	pressures = 'linear fit: gradient {:E}  intercept {:E}, atomic units'.format(fit[0],fit[1])
	data.append(pressures)
	#Convert eV to Joules
	eV_to_Joules = 1.60217646E-19
	#Convert Angstrom to cubic metre
	angstrom_to_metre = 1E-10
	x = x*(angstrom_to_metre**3)
	y = y*eV_to_Joules
	fit = numpy.polyfit(x, y, 1)
	pressures = 'linear fit: gradient {:E}  intercept {:E}, SI units'.format(fit[0],fit[1])
	data.append(pressures)
	EVQ_pressure_correction = 'Pressure correction: {:E} GPa'.format(fit[0]*1E-9)
	data.append(EVQ_pressure_correction)
	return data

def output_ZnOX_positions(atoms):
	#returns the distances from each atom to the atom number specified
	data = []
	heading = 'ZnOX_positions_output'
	input = os.path.realpath('OUTCAR')
	data.append(heading)
	data.append(input)
	heading = 'ZnOX_cell_parameters'
	data.append(heading)
	columns = '{:<18}{:^18}{:^18}{:^18}{:^18}{:^18}{:^18}'.format('Volume','a','b','c','alpha','beta','gamma')
	data.append(columns)
	asetools.cellparam(atoms)
	cell_data = '{:<18}{:^18}{:^18}{:^18}{:^18}{:^18}{:^18}'.format(atoms.get_volume(),asetools.cellparam(atoms)[0],asetools.cellparam(atoms)[1],asetools.cellparam(atoms)[2],asetools.cellparam(atoms)[3],asetools.cellparam(atoms)[4],asetools.cellparam(atoms)[5])
	data.append(cell_data)
	heading = 'ZnOX_cell_matrix'
	data.append(heading)
	columns = '{:<18}{:^18}{:^18}'.format(atoms.get_cell()[0][0],atoms.get_cell()[0][1],atoms.get_cell()[0][2])
	data.append(columns)
	columns = '{:<18}{:^18}{:^18}'.format(atoms.get_cell()[1][0],atoms.get_cell()[1][1],atoms.get_cell()[1][2])
	data.append(columns)
	columns = '{:<18}{:^18}{:^18}'.format(atoms.get_cell()[2][0],atoms.get_cell()[2][1],atoms.get_cell()[2][2])
	data.append(columns)
	heading = 'XO4_tetrahedra_positions_and_bondlengths'
	data.append(heading)
	columns = '{:<8}{:^15}{:^18}{:^18}{:^18}{:^18}'.format('Atom','Atom_number','O-X_bond_length','Atom_a_coord','Atom_b_coord','Atom_c_coord')
	data.append(columns)
	ions = [71,42,61,59,53]
	# loop over ions
	for i in range(len(ions)):
		n = ions[i]
		atom = atoms.get_chemical_symbols()[n]
		bond_length = atoms.get_distance(71,n,mic=1)
		coord_a = atoms.get_positions()[n,0]
		coord_b = atoms.get_positions()[n,1]
		coord_c = atoms.get_positions()[n,2]
		info = '{:<8}{:^15}{:^18}{:^18}{:^18}{:^18}'.format(atom,n+1,bond_length,coord_a,coord_b,coord_c)
		data.append(info)
	outfile = 'ZnOX_positions.txt'
	fout = open(outfile,'w')
	for item in data:
		print>>fout, item
	fout.close
	return(data)

def output_ZnOX_X_tetrahedra(atoms):
	"""returns the positions and a POSCAR showing the tetrahedra around the X defect in the ZnOX project
	specify an ASE atoms object to work on and it has to be setup as the ZnOX positions were done,
	it gets the tetrahedra made by the 72nd ion (X) and surrounding Ox atoms - 42, 61, 59 and 53
	returns a POSCAR formatted file with only the defective tetrahedron
	i.e. XO4^-4 tetrahedron in the supercell"""
	del atoms[[list(range(0,42)) + list(range(43,53)) + list(range(54,59)) + list(range(60,61)) + list(range(62,71))]]
	filename = 'POSCAR_ZnOX_' + atoms.get_chemical_symbols()[4] + '_defect_tetradron'
	ase.io.write(filename, atoms, format='vasp')

def output_distances_to_point(atoms,a,b,c):
	"""returns the distances from each atom to the fractional coordinate specified"""
	data = []
	heading = 'distances_to_point'
	subheading = 'origin'
	origin_col = '{:<10} {:^10} {:^10}'.format('a','b','c')
	origin = '{:<10} {:^10} {:^10}'.format(a,b,c)
	input = os.path.realpath('OUTCAR')
	columns = '{:<8} {:^15} {:^18} {:^18} {:^18} {:^18}'.format('Atom','Atom_number','O-X_bond_length','Atom_a_coord','Atom_b_coord','Atom_c_coord')
	data.append(heading)
	data.append(input)
	data.append(subheading)
	data.append(origin_col)
	data.append(origin)
	data.append(columns)
	# Add a Xe atom in the coordinate to measure from
	new = ase.Atoms('Xe')
	new.set_cell(atoms.get_cell())
	new.set_scaled_positions([a,b,c])
	atoms.append('Xe')
	atoms[(len(atoms.get_positions())-1)].set_position([new.get_positions()[0][0],new.get_positions()[0][1],new.get_positions()[0][2]])
	# Create a table of distances
	point = len(atoms.get_chemical_symbols())-1
	for n in range(len(atoms.get_chemical_symbols())-1):
		atom = atoms.get_chemical_symbols()[n]
		bond_length = atoms.get_distance(point,n,mic=1)
		coord_a = atoms.get_positions()[n,0]
		coord_b = atoms.get_positions()[n,1]
		coord_c = atoms.get_positions()[n,2]
		info = '{:<8} {:^15} {:^18} {:^18} {:^18} {:^18}'.format(atom,n+1,bond_length,coord_a,coord_b,coord_c)
		data.append(info)
	outfile = 'ZnOX_distances_to_point_{}-{}-{}.txt'.format(a,b,c)
	fout = open(outfile,'w')
	for item in data:
		print>>fout, item
	fout.close
	del atoms[len(atoms.get_chemical_symbols())-1]
	return(data)

def output_distances_to_atom(atoms,atom):
	"""returns the distances and position of the nearest n atoms to the fractional coordinate specified"""
	data = []
	heading = 'distances_to_atom'
	subheading = 'origin'
	a = atoms.get_positions()[atom-1][0]
	b = atoms.get_positions()[atom-1][1]
	c = atoms.get_positions()[atom-1][2]
	origin_col = '{:<10} {:^10} {:^10}'.format('a','b','c')
	origin = '{:<10} {:^10} {:^10}'.format(a,b,c)
	input = os.path.realpath('OUTCAR')
	columns = '{:<8} {:^15} {:^18} {:^18} {:^18} {:^18}'.format('Atom','Atom_number','O-X_bond_length','Atom_a_coord','Atom_b_coord','Atom_c_coord')
	data.append(heading)
	data.append(input)
	data.append(subheading)
	data.append(origin_col)
	data.append(origin)
	data.append(columns)
	# Create a table of distances
	point = atom-1
	for n in range(len(atoms.get_chemical_symbols())):
		species = atoms.get_chemical_symbols()[n]
		bond_length = atoms.get_distance(point,n,mic=1)
		coord_a = atoms.get_positions()[n,0]
		coord_b = atoms.get_positions()[n,1]
		coord_c = atoms.get_positions()[n,2]
		info = '{:<8} {:^15} {:^18} {:^18} {:^18} {:^18}'.format(species,n+1,bond_length,coord_a,coord_b,coord_c)
		data.append(info)
	outfile = 'ZnOX_distances_to_atom{number}_{symbol}.txt'.format(number=atom,symbol=atoms.get_chemical_symbols()[atom-1])
	fout = open(outfile,'w')
	for item in data:
		print>>fout, item
	fout.close
	return(data)

def output_magnetic_moments_vs_distance_to_atom(atoms,atom,outcar_file):
	"""returns the distances, magnetic moment and position of the nearest n atoms to the fractional coordinate specified
		give it an "atoms" object, the atom number, and an 'OUTCAR' file"""
	import re
	from StringIO import StringIO
	data = []
	heading = 'magnetic_moments_and_distances_to_atom'
	subheading = 'origin'
	a = atoms.get_positions()[atom-1][0]
	b = atoms.get_positions()[atom-1][1]
	c = atoms.get_positions()[atom-1][2]
	origin_col = '{:<18} {:^18} {:^18}'.format('a','b','c')
	origin = '{:<18} {:^18} {:^18}'.format(a,b,c)
	input = os.path.realpath('OUTCAR')
	columns = '{:<6} {:^13} {:^18} {:^18} {:^18} {:^18} {:^18}'.format('Atom','Atom_number','Atom-X_distance','Magnetic_moment','Atom_a_coord','Atom_b_coord','Atom_c_coord')
	data.append(heading)
	data.append(input)
	data.append(subheading)
	data.append(origin_col)
	data.append(origin)
	data.append(columns)
	f = open(outcar_file)
	magnetisation_data = []
	for line in f:
		if "magnetization (x)" in line:
			for i in range(len(atoms.get_positions())+4):
				item = next(f)
				magnetisation_data.append(item)
	f.close
#	Create a table of distances
	point = atom-1
	for n in range(len(atoms.get_chemical_symbols())):
		species = atoms.get_chemical_symbols()[n]
		bond_length = atoms.get_distance(point,n,mic=1)
		coord_a = atoms.get_positions()[n,0]
		coord_b = atoms.get_positions()[n,1]
		coord_c = atoms.get_positions()[n,2]
		magnetisation = magnetisation_data[n+3].split()
		info = '{:<6} {:^13} {:^18} {:^18} {:^18} {:^18} {:^18}'.format(species,n+1,bond_length,magnetisation[-1],coord_a,coord_b,coord_c)
		data.append(info)
	outfile = 'Magnetic_moments_on_atoms_surrounding_atom{number}_{symbol}.txt'.format(number=atom,symbol=atoms.get_chemical_symbols()[atom-1])
	fout = open(outfile,'w')
	for item in data:
		print>>fout, item
	fout.close

def defect_strain_matrices(label,conc,Latt):
	from ase.calculators.vasp import Vasp
	calc = Vasp(restart=1)
	atoms = calc.get_atoms()
	LattDef = numpy.array(atoms.get_cell())
	# LattDef = correct_lattice_matrix(numpy.array(atoms.get_cell()))
	dstrain = (1/conc)*numpy.dot(numpy.subtract(LattDef,Latt),numpy.linalg.inv(Latt))
	dstrain_vec = asetools.tensorise(dstrain)
	outfile = 'defect_strain.txt'
	fout = open(outfile,'w')
	print>>fout, label
	print>>fout, "concentration {:^10}".format(conc)
	print>>fout, "defect strain matrix"
	for item in dstrain:
		print>>fout, '   '.join(map(str, item))
	print>>fout, "defect strain vector"
	for item in dstrain_vec:
		print>>fout, item
	print>>fout, "Defect supercell lattice"
	for item in LattDef:
		print>>fout, '   '.join(map(str, item))
	print>>fout, "Perfect supercell lattice"
	for item in Latt:
		print>>fout, '   '.join(map(str, item))
	fout.close

def correct_lattice_matrix(matrix):
	"""adjusts the lattice matrix to effectively flip the a and b cell parameters"""
	correct_matrix = numpy.copy(matrix)
	correct_matrix[0,0] = matrix[1,1]
	correct_matrix[0,1] = matrix[1,0]
	correct_matrix[1,0] = matrix[0,1]
	correct_matrix[1,1] = matrix[0,0]
	return correct_matrix

def output_cell_params():
	import fnmatch
	from ase.calculators.vasp import Vasp
	outfile = 'cell_params.txt'
	fout = open(outfile,'w')
	prelim_matches = []
	matches = []
	for root, dirnames, filenames in os.walk('.'):
		for filename in fnmatch.filter(filenames, 'OUTCAR'):
			prelim_matches.append(os.path.abspath(root))
	for dir in prelim_matches:
		if os.path.isfile(os.path.join(dir,"CONTCAR")):
			matches.append(dir)
	for dir in matches:
		os.chdir(dir)
		label = os.path.abspath(".")
		print>>fout, label
		print>>fout, "{:<14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} ".format('volume','a','b','c','alpha','beta','gamma')
		calc = Vasp(restart=1)
		atoms = calc.get_atoms()
		params = asetools.cellparam(atoms)
		volume = atoms.get_volume()
		print>>fout, "{:<14} {:^14} {:^14} {:^14} {:^14} {:^14} {:^14} ".format(volume,params[0],params[1],params[2],params[3],params[4],params[5])
	fout.close

def get_potential_correction_ZnOX(data_file,line_number):
	fp = open(data_file)
	for i, line in enumerate(fp):
		if i == (line_number-1):
			numbers = line.split()
		elif i > line_number:
			break
	fp.close()
	return numbers[2]

def edit_vasprun_for_DOS(file,nkpts):
	""" edits the VASP xml file to let P4vasp plot it directly as a DOS plot """
	from lxml import etree
	kpoints_to_remove = nkpts
	file = open('vasprun.xml', "r")
	elem = etree.parse(file)
	#remove excess kpoints
	for N in range(1,(kpoints_to_remove+1)):
		for i in elem.xpath('/modeling/calculation/eigenvalues/array/set/set[@comment="spin 1"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/eigenvalues/array/set/set[@comment="spin 2"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/projected/eigenvalues/array/set/set[@comment="spin 1"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/projected/eigenvalues/array/set/set[@comment="spin 2"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/projected/array/set/set[@comment="spin1"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/projected/array/set/set[@comment="spin2"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/projected/array/set/set[@comment="spin 1"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)
		for i in elem.xpath('/modeling/calculation/projected/array/set/set[@comment="spin 2"]/set[@comment="kpoint %i"]' % N) :
			i.getparent().remove(i)

	# reset index numbers of kpoints
	for comment in elem.xpath('/modeling/calculation/eigenvalues/array/set/*/set'):
		old_string = (comment.attrib['comment']).split()
		number = int(old_string[1])
		comment.attrib['comment'] = 'kpoint %i' % (number - nkpts)
	for comment in elem.xpath('/modeling/calculation/projected/eigenvalues/array/set/*/set'):
		old_string = (comment.attrib['comment']).split()
		number = int(old_string[1])
		comment.attrib['comment'] = 'kpoint %i' % (number - nkpts)
	for comment in elem.xpath('/modeling/calculation/projected/array/set/*/set'):
		old_string = (comment.attrib['comment']).split()
		number = int(old_string[1])
		comment.attrib['comment'] = 'kpoint %i' % (number - nkpts)
	file.close

	outfile = 'DOS_vasprun.xml'
	fout = open(outfile,'w')
	elem.write(fout)
	fout.close

	file = open(outfile,'r')
	kpointlist = 'varray name="kpointlist"'
	weights = 'varray name="weights"'
	for i, line in enumerate(file, 1):
		if kpointlist in line:
			kp_start = i
		if weights in line:
			wt_start = i
	file.close
	print "Now run this sed command from the shell:"
	print "sed -i -e '{0},{1}{2}' -e '{3},{4}{5}' {6}".format(wt_start+1, wt_start+nkpts, "d", kp_start+1, kp_start+nkpts, "d", os.path.abspath(outfile))




