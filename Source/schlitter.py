#!/usr/local/bin/python

"""This is the main file to compute
the intramolecular entropy based on
Schlitter's method
"""

import os
import sys
import gzip
import argparse
import numpy as np
from utilities import set_units, kabsch, form_disp_matrix
from ParseLammpsDataFile import read_preliminary_data, read_atomic_masses, read_atomic_info, read_bonds
from ParseLammpsDumpFile import read_configuration

def intramolecular_entropy():
 """This function is coordinating the calculation of the intramolecular entropy """

 parser = argparse.ArgumentParser(description='A python script to compute the Intramolecular Entropy based on the Schlitter\'s method. Entropies are computed in Boltzmann\'s constant units, i.e. they are normalized by k_B. The system is assumed to be periodic in all three directions.')

 parser.add_argument("Temperature", type=float, help="The temperature of the system in Kelvin")
 parser.add_argument("OutputEntropyFile", help="The name of the output file where the eigenvalues of the covariance matrix will be stored as well as their contribution to the intramolecular entropy. The name should include the full path to the file. The file is assumed to be new.")
 parser.add_argument("LAMMPSInputDataFile", help="The name of the lammps input data file. The atom_style is assumed to be \"full\" otherwise the information about molecules cannot be retrieved from the data file.")
 parser.add_argument("DumpFileType", help="Define the type of the dump file. It can be either \"xyz\" which means that the file is an xyz one, 'scaled' which means the coordinates are scaled and the dump is the usual lammps file or 'typical' which is a regular LAMMPS dump file with the format \"id x y z\".")
 parser.add_argument("LAMMPSInputDumpFile", nargs='+', help="A list with the name of the lammps input trajectory files. They should include a full path.")
 parser.add_argument("--units", help="define the units to be used in the analysis. The available unit systems are SI, CGS, METAL, LJ, MICRO, NANO, REAL and ELECTRON. It should be written in capital letters. If the specified unit systems is not implemented then the analysis will be stopped immediately and an error message will be printed in the screen. The default unit system is REAL.")
 parser.add_argument("--restart", help="If specified it means that the analysis will restart after the unfolding.")

 args = parser.parse_args()

 # open an intermediate dump file to store the unfolded coordinates
 unfolded_file = open("UnfoldedDumpFile.txt", "ab+")

 # open an intermediate dump file to store the unfolded coordinates with removed rotational & translational effects
 intermediate_file = open("IntermediateDumpFile.txt", "ab+")

 if args.units is None:
        reduced_planck, boltzmann_constant = set_units("REAL")
 else:
        reduced_planck, boltzmann_constant = set_units(args.units)

 beta = 1 / (boltzmann_constant*args.Temperature)

 # open the LAMMPS data file to read the molecules
 with open(args.LAMMPSInputDataFile, "r") as lammps_data_file:
        num_atoms, num_bonds, num_atom_types, xboxlength, yboxlength, zboxlength = read_preliminary_data(lammps_data_file)

        atom_coord = np.zeros((num_atoms, 3))
        atom_type = np.zeros(num_atoms, dtype=int)
        dispersity, molecule_id = [0] * num_atoms, [0] * num_atoms
        table_mass, atom_mass = np.zeros(num_atom_types), np.zeros(num_atoms)
        ref_atom = np.zeros((num_atoms, 3))
        mean_atom_pos = np.zeros((num_atoms, 3))
        sorted_bond_sequence = []

        read_atomic_masses(lammps_data_file, table_mass, num_atom_types)
        read_atomic_info(lammps_data_file, num_atoms, molecule_id, atom_type, dispersity)
        read_bonds(lammps_data_file, num_bonds, sorted_bond_sequence)

 np_mol_id = np.array(molecule_id, dtype=int)

 num_molecules = max(molecule_id) + 1
 dispersity = dispersity[:num_molecules]
 molecule_mass, min_rg = np.zeros(num_molecules), np.zeros(num_molecules)

 # find the maximum number of atoms that one molecule contains
 max_atoms_per_mol = max(dispersity)
 disp_matrix = np.zeros(shape=(num_molecules, 3*max_atoms_per_mol, 3*max_atoms_per_mol))

 for itype in range(0, num_atom_types):
    atom_mass[atom_type == itype] = table_mass[itype]

 sq_mass = np.sqrt(atom_mass)

 atoms_in_mol = []
 for imol in range(0, num_molecules): # compute the COM  and the squared radius of gyration of each molecule
  atoms_in_mol.append(np.where(np_mol_id == imol)[0])
  molecule_mass[imol] = np.sum(atom_mass[np_mol_id == imol])

 unfold_left, unfold_right = [], []
 for imol in range(0, num_molecules):
  exclude, current = [], []
  k = molecule_id.index(imol)
  current.append(k) # the starting atom is appended to the current list

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

 num_confs = 0
 if args.restart is None or args.restart < 1:

  for dump_file in args.LAMMPSInputDumpFile:

   cur_dump_file = gzip.open(dump_file, "rb")

   while True: # perform the calculations by analyzing all available configurations

    if read_configuration(args.DumpFileType, cur_dump_file, num_atoms, atom_coord):
        break

    num_confs += 1

    for iat, jat in zip(unfold_left, unfold_right):
      atom_coord[jat, 0] -= xboxlength*round((atom_coord[jat, 0] - atom_coord[iat, 0])/xboxlength)
      atom_coord[jat, 1] -= yboxlength*round((atom_coord[jat, 1] - atom_coord[iat, 1])/yboxlength)
      atom_coord[jat, 2] -= zboxlength*round((atom_coord[jat, 2] - atom_coord[iat, 2])/zboxlength)

    for imol in range(0, num_molecules):
     com_x = np.sum(np.multiply(atom_mass[atoms_in_mol[imol]], atom_coord[atoms_in_mol[imol], 0])) / molecule_mass[imol]
     com_y = np.sum(np.multiply(atom_mass[atoms_in_mol[imol]], atom_coord[atoms_in_mol[imol], 1])) / molecule_mass[imol]
     com_z = np.sum(np.multiply(atom_mass[atoms_in_mol[imol]], atom_coord[atoms_in_mol[imol], 2])) / molecule_mass[imol]

     atom_coord[atoms_in_mol[imol], 0] -= com_x
     atom_coord[atoms_in_mol[imol], 1] -= com_y
     atom_coord[atoms_in_mol[imol], 2] -= com_z

     square_rg = np.sum(atom_mass[atoms_in_mol[imol]]*(atom_coord[atoms_in_mol[imol], 0]**2 + \
      atom_coord[atoms_in_mol[imol], 1]**2 + atom_coord[atoms_in_mol[imol], 2]**2)) / molecule_mass[imol]

     if num_confs == 1 or np.sqrt(square_rg) < min_rg[imol]:
                        ref_atom[atoms_in_mol[imol], 0] = atom_coord[atoms_in_mol[imol], 0].copy()
                        ref_atom[atoms_in_mol[imol], 1] = atom_coord[atoms_in_mol[imol], 1].copy()
                        ref_atom[atoms_in_mol[imol], 2] = atom_coord[atoms_in_mol[imol], 2].copy()
                        min_rg[imol] = np.sqrt(square_rg)

    np.savetxt(unfolded_file, np.c_[atom_coord[:, 0], atom_coord[:, 1], atom_coord[:, 2]])

  with open("ReferenceConfiguration.txt", "wb") as ref_state_file:
            for imol in range(0, num_molecules):
                np.savetxt(ref_state_file, np.c_[atoms_in_mol[imol], \
                ref_atom[atoms_in_mol[imol], 0], ref_atom[atoms_in_mol[imol], 1], ref_atom[atoms_in_mol[imol]], 2])

 else:
        with open("ReferenceConfiguration.txt", "rb") as ref_state_file:
            atom_id, ref_atom[:, 0], ref_atom[:, 1], ref_atom[:, 2] = \
                   np.loadtxt(ref_state_file, dtype='int, float, float, float', unpack=True)
            ref_atom[:, 0] = ref_atom[atom_id, 0]
            ref_atom[:, 1] = ref_atom[atom_id, 1]
            ref_atom[:, 2] = ref_atom[atom_id, 2]

 unfolded_file.seek(0) # rewind the intermediate dump file

 if num_confs < (3*max_atoms_per_mol + 1):
        sys.exit('The number of configurations is smaller than the maximum number of atoms in the system')

 iconf = 0
 while True:

        if read_configuration('intermediate', unfolded_file, num_atoms, atom_coord):
            break

        iconf += 1
        for imol in range(0, num_molecules):
            atom_coord[atoms_in_mol[imol], :] = \
                        kabsch(atom_coord[atoms_in_mol[imol], :], ref_atom[atoms_in_mol[imol], :])

        np.savetxt(intermediate_file, np.c_[atom_coord[:, 0], atom_coord[:, 1], atom_coord[:, 2]])

        mean_atom_pos[:, 0] += (atom_coord[:, 0] - mean_atom_pos[:, 0]) / float(iconf)
        mean_atom_pos[:, 1] += (atom_coord[:, 1] - mean_atom_pos[:, 1]) / float(iconf)
        mean_atom_pos[:, 2] += (atom_coord[:, 2] - mean_atom_pos[:, 2]) / float(iconf)

 unfolded_file.close()
 os.remove("UnfoldedDumpFile.txt") # delete the intermediate file with the unfolded coordinates

 intermediate_file.seek(0) # rewind the intermediate dump file

 iconf = 0
 while True: # perform the calculations of the displacement matrix by analyzing all available configurations

        if read_configuration('intermediate', intermediate_file, num_atoms, atom_coord):
            break

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

 intermediate_file.close()
 os.remove("IntermediateDumpFile.txt")

 with open(args.OutputEntropyFile, "wb") as output_file:
        for imol in range(0, num_molecules):
            eigval = np.linalg.eigvalsh(disp_matrix[imol][0:3*dispersity[imol]][0:3*dispersity[imol]])
            sort_eigenvalues = eigval.argsort()
            omega = 1 / np.sqrt(eigval[sort_eigenvalues[6: ]]*beta)
            mode_entropy = (beta*omega*reduced_planck)/(np.exp(beta*omega*reduced_planck)-1) \
                                                   - np.log(1 - np.exp(-beta*omega*reduced_planck))
            output_file.write("The total intramolecular entropy is {} \n".format(sum(mode_entropy)))
            np.savetxt(output_file, np.c_[eigval[sort_eigenvalues[6: ]], mode_entropy])

if __name__ == "__main__":
    intramolecular_entropy()
