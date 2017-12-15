#!/usr/local/bin/python

"""This is the main file to compute
the intramolecular entropy based on
Schlitter's method
"""

from os import remove
from sys import exit
from gzip import open as Open
from argparse import ArgumentParser
from numpy import zeros, sqrt, amax, count_nonzero, where, multiply, savetxt, c_, array, loadtxt, linalg, exp, log, sum as Sum
from utilities import set_units, kabsch, form_disp_matrix
from parse_data_file import read_preliminary_data, read_atomic_masses, read_atomic_info, read_bonds
from parse_dump_file import read_configuration

def intramolecular_entropy():
 """This function is coordinating the calculation of the intramolecular entropy """

 parser = ArgumentParser(description='A python script to compute the Intramolecular Entropy based on the Schlitter\'s method. Entropies are computed in Boltzmann\'s constant units, i.e. they are normalized by k_B. The system is assumed to be periodic in all three directions.')

 parser.add_argument("Temperature", type=float, help="The temperature of the system in Kelvin")
 parser.add_argument("OutputEntropyFile", help="The name of the output file where the eigenvalues of the covariance matrix will be stored as well as their contribution to the intramolecular entropy. The name should include the full path to the file. The file is assumed to be new.")
 parser.add_argument("InputDataFile", help="The name of the lammps input data file. The atom_style is assumed to be \"full\" otherwise the information about molecules cannot be retrieved from the data file.")
 parser.add_argument("DumpFileType", help="Define the type of the dump file. It can be either \"xyz\" which means that the file is an xyz one, 'scaled' which means the coordinates are scaled and the dump is the usual lammps file or 'typical' which is a regular LAMMPS dump file with the format \"id x y z\".")
 parser.add_argument("InputDumpFile", nargs='+', help="A list with the name of the lammps input trajectory files. They should include a full path.")
 parser.add_argument("--units", help="define the units to be used in the analysis. The available unit systems are SI, CGS, METAL, LJ, MICRO, NANO, REAL and ELECTRON. It should be written in capital letters. If the specified unit systems is not implemented then the analysis will be stopped immediately and an error message will be printed in the screen. The default unit system is REAL.")
 parser.add_argument("--pbc", nargs=3, default=[True, True, True], help="define the periodicity of the simulation box. Three boolean logical values True or False must be provided for the x-, y- and z- axis.")
 parser.add_argument("--restart", help="If specified it means that the analysis will restart after the unfolding.")
 parser.add_argument("--yasp", help="If specified it means that the input files have been created by YASP code")

 args = parser.parse_args()

 if args.yasp is True:
        args.units = "YASP"
        args.DumpFileType = "xyz"

 if args.units is None:
        reduced_planck, boltzmann_constant = set_units("REAL")
 else:
        reduced_planck, boltzmann_constant = set_units(args.units)

 beta = 1 / (boltzmann_constant*args.Temperature)

 # Open the data file to read the molecules
 sorted_bond_sequence = []
 with Open(args.InputDataFile, "r") as data_file:
        num_atoms, num_bonds, num_atom_types, xboxlength, yboxlength, zboxlength = read_preliminary_data(data_file)

        atom_type = zeros(num_atoms, dtype=int)
        molecule_id = zeros(num_atoms, dtype=int)
        atom_mass = zeros(num_atoms)

        read_atomic_info(data_file, num_atoms, molecule_id, atom_type)
        read_atomic_masses(data_file, atom_mass, atom_type, num_atom_types)

        read_bonds(data_file, num_bonds, sorted_bond_sequence)

 atom_coord = zeros((num_atoms, 3))
 sq_mass = sqrt(atom_mass)
 ref_atom = zeros((num_atoms, 3))
 mean_atom_pos = zeros((num_atoms, 3))

 num_molecules = amax(molecule_id) + 1
 dispersity = [0] * num_molecules
 for imol in range(0, num_molecules):
     dispersity[imol] = count_nonzero(molecule_id == imol)
 molecule_mass, min_rg = zeros(num_molecules), zeros(num_molecules)

 # find the maximum number of atoms that one molecule contains
 max_atoms_per_mol = max(dispersity)
 disp_matrix = zeros(shape=(num_molecules, 3*max_atoms_per_mol, 3*max_atoms_per_mol))

 atoms_in_mol = []
 for imol in range(0, num_molecules):
  atoms_in_mol.append(where(molecule_id == imol)[0])
  molecule_mass[imol] = Sum(atom_mass[molecule_id == imol])

 unfold_left, unfold_right = [], []
 for imol in range(0, num_molecules):
  exclude, current = [], []
  current.append(atoms_in_mol[imol][0]) # the starting atom is appended to the current list

  while len(exclude) < dispersity[imol]:
   new_list = []
   for jat in current:
      if jat not in exclude:
         exclude.append(jat)

      for lat in sorted_bond_sequence[jat]:
         if lat not in exclude:
            new_list.append(lat)
            unfold_left.append(jat)
            unfold_right.append(lat)
   # move the contents of the new list to the current list
   current = new_list

 with Open("UnfoldedDumpFile.txt.gz", "ab+") as unfolded_file:
  num_confs = 0
  if args.restart is None:

   for dump_file in args.InputDumpFile:

    with Open(dump_file, "rb") as cur_dump_file:

    # perform the calculations by analyzing all available configurations
     while not read_configuration(args.DumpFileType, cur_dump_file, num_atoms, atom_coord):

      num_confs += 1
      print("The current timestep is {}".format(num_confs))

      for iat, jat in zip(unfold_left, unfold_right):
       atom_coord[jat, 0] -= xboxlength*round((atom_coord[jat, 0] - atom_coord[iat, 0])/xboxlength)
       atom_coord[jat, 1] -= yboxlength*round((atom_coord[jat, 1] - atom_coord[iat, 1])/yboxlength)
       atom_coord[jat, 2] -= zboxlength*round((atom_coord[jat, 2] - atom_coord[iat, 2])/zboxlength)

      for imol in range(0, num_molecules):
       com_x = Sum(multiply(atom_mass[atoms_in_mol[imol]], \
               atom_coord[atoms_in_mol[imol], 0])) / molecule_mass[imol]
       com_y = Sum(multiply(atom_mass[atoms_in_mol[imol]], \
               atom_coord[atoms_in_mol[imol], 1])) / molecule_mass[imol]
       com_z = Sum(multiply(atom_mass[atoms_in_mol[imol]], \
               atom_coord[atoms_in_mol[imol], 2])) / molecule_mass[imol]

       atom_coord[atoms_in_mol[imol], 0] -= com_x
       atom_coord[atoms_in_mol[imol], 1] -= com_y
       atom_coord[atoms_in_mol[imol], 2] -= com_z

       square_rg = Sum(atom_mass[atoms_in_mol[imol]]*(atom_coord[atoms_in_mol[imol], 0]**2 + \
        atom_coord[atoms_in_mol[imol], 1]**2 + atom_coord[atoms_in_mol[imol], 2]**2)) / molecule_mass[imol]

       if num_confs == 1 or sqrt(square_rg) < min_rg[imol]:
                        ref_atom[atoms_in_mol[imol], 0] = atom_coord[atoms_in_mol[imol], 0].copy()
                        ref_atom[atoms_in_mol[imol], 1] = atom_coord[atoms_in_mol[imol], 1].copy()
                        ref_atom[atoms_in_mol[imol], 2] = atom_coord[atoms_in_mol[imol], 2].copy()
                        min_rg[imol] = sqrt(square_rg)

      savetxt(unfolded_file, c_[atom_coord[:, 0], atom_coord[:, 1], atom_coord[:, 2]], fmt='%.10e')

    #remove(dump_file)

   if num_confs < (3*max_atoms_per_mol + 1):
        exit('The number of configurations is smaller than the maximum number of atoms in the system')


   with Open("ReferenceConfiguration.txt", "wb") as ref_state_file:
            for imol in range(0, num_molecules):
                savetxt(ref_state_file, c_[array(atoms_in_mol[imol]), \
                ref_atom[array(atoms_in_mol[imol]), 0], ref_atom[array(atoms_in_mol[imol]), 1], ref_atom[array(atoms_in_mol[imol]), 2]])

  else:
        with Open("ReferenceConfiguration.txt", "rb") as ref_state_file:
            atom_id, ref_atom[:, 0], ref_atom[:, 1], ref_atom[:, 2] = \
                   loadtxt(ref_state_file, dtype='float, float, float, float', unpack=True)
            atom_id = atom_id.astype(int)
            ref_atom[:, 0] = ref_atom[atom_id, 0]
            ref_atom[:, 1] = ref_atom[atom_id, 1]
            ref_atom[:, 2] = ref_atom[atom_id, 2]


 with Open("IntermediateDumpFile.txt.gz", "wb+") as intermediate_file:
  with Open("UnfoldedDumpFile.txt.gz", "rb+") as unfolded_file:

   iconf = 0
   while not read_configuration('intermediate', unfolded_file, num_atoms, atom_coord):

        print("step {}".format(iconf))
        iconf += 1
        for imol in range(0, num_molecules):
            atom_coord[atoms_in_mol[imol], :] = \
                        kabsch(atom_coord[atoms_in_mol[imol], :], ref_atom[atoms_in_mol[imol], :])

        savetxt(intermediate_file, \
                         c_[atom_coord[:, 0], atom_coord[:, 1], atom_coord[:, 2]], fmt='%.10e')

        mean_atom_pos[:, 0] += (atom_coord[:, 0] - mean_atom_pos[:, 0]) / float(iconf)
        mean_atom_pos[:, 1] += (atom_coord[:, 1] - mean_atom_pos[:, 1]) / float(iconf)
        mean_atom_pos[:, 2] += (atom_coord[:, 2] - mean_atom_pos[:, 2]) / float(iconf)

 remove("UnfoldedDumpFile.txt.gz") # delete the intermediate file with the unfolded coordinates

 with Open("IntermediateDumpFile.txt.gz", "rb+") as intermediate_file:

  iconf = 0
  while not read_configuration('intermediate', intermediate_file, num_atoms, atom_coord):

        print("The current timestep is {}".format(iconf))

        atom_coord[:, 0] -= mean_atom_pos[:, 0]
        atom_coord[:, 1] -= mean_atom_pos[:, 1]
        atom_coord[:, 2] -= mean_atom_pos[:, 2]

        iconf += 1
        # form the mass-weighted displacement matrix for each molecule
        for imol in range(0, num_molecules):
            nat = dispersity[imol]
            displacement_matrix = form_disp_matrix(atom_coord[atoms_in_mol[imol], :], sq_mass[atoms_in_mol[imol]])
            disp_matrix[imol][0:3*nat][0:3*nat] += \
             (displacement_matrix - disp_matrix[imol][0:3*nat][0:3*nat])/ float(iconf)

 remove("IntermediateDumpFile.txt.gz")

 with Open(args.OutputEntropyFile, "w") as output_file:
        for imol in range(0, num_molecules):
            eigval = linalg.eigvalsh(disp_matrix[imol][0:3*dispersity[imol]][0:3*dispersity[imol]])
            sort_eigenvalues = eigval.argsort()
            omega = 1 / sqrt(eigval[sort_eigenvalues[6: ]]*beta)
            mode_entropy = (beta*omega*reduced_planck)/(exp(beta*omega*reduced_planck)-1) \
                                                   - log(1 - exp(-beta*omega*reduced_planck))
            output_file.write("The total intramolecular entropy is {} \n".format(sum(mode_entropy)))
            for ival, jval in zip(eigval[sort_eigenvalues[6: ]], mode_entropy):
              output_file.write("{} {} \n".format(ival, jval))

if __name__ == "__main__":
    intramolecular_entropy()
