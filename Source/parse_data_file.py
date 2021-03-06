#!/usr/local/bin/python

"""
This module process the input LAMMPS data file.
It contains the functions:
-read_preliminary_data
-read_atomic_masses
-read_atomic_info
-read_bonds
"""

from numpy import zeros, amax, concatenate

def read_preliminary_data(input_file):
    """ find how many atoms & atom types are contained in the system"""

    sim_box = zeros((3, 3))

    for line in input_file:
        temp = line.split()
        if len(temp) > 1:
            if temp[1] == 'atoms':
                num_atoms = int(temp[0])
            elif temp[1] == 'bonds':
                num_bonds = int(temp[0])
            elif len(temp) > 2:
                if temp[1] == 'atom' and temp[2] == 'types':
                    num_atom_types = int(temp[0])
                elif len(temp) > 3:
                    if temp[2] == "xlo" and temp[3] == "xhi":
                        sim_box[0, 0] = float(temp[1]) - float(temp[0])
                    elif temp[2] == "ylo" and temp[3] == "yhi":
                        sim_box[1, 1] = float(temp[1]) - float(temp[0])
                    elif temp[2] == "zlo" and temp[3] == "zhi":
                        sim_box[2, 2] = float(temp[1]) - float(temp[0])
                    elif len(temp) > 5:
                        if temp[3] == "xy" and temp[4] == "xz" and temp[5] == "yz":
                            sim_box[0, 1] = float(temp[0])
                            sim_box[0, 2] = float(temp[1])
                            sim_box[1, 2] = float(temp[2])

    input_file.seek(0)  # rewind the LAMMPS data file

    return num_atoms, num_bonds, num_atom_types, sim_box

def read_atomic_masses(input_file, atom_mass, atom_type, number_atom_types):
    """ find the masses of all atom types in the system"""

    table_mass = zeros(number_atom_types)

    while True:
        line = input_file.readline()
        temp = line.split()
        if temp:
            if temp[0] == 'Masses':
                while True:
                    line = input_file.readline()
                    temp = line.split()
                    if temp:
                        table_mass[int(temp[0])-1] = float(temp[1])
                        for _ in range(1, number_atom_types):
                            line = input_file.readline()
                            temp = line.split()
                            table_mass[int(temp[0])-1] = float(temp[1])

                        # assign the mass to each atom here
                        for itype in range(0, number_atom_types):
                            atom_mass[atom_type == itype] = table_mass[itype]

                        # rewind the LAMMPS data file
                        input_file.seek(0)
                        return  

def read_atomic_info(input_file, num_atoms, molecule_id, atom_type):
    """It extracts the molecule id and atom type of the input structure """
    while True:
        line = input_file.readline()
        temp = line.split()
        if temp:
            if temp[0] == 'Atoms':
                while True:
                    line = input_file.readline()
                    temp = line.split()
                    if temp:
                        iat = int(temp[0]) - 1
                        imol = int(temp[1]) - 1
                        molecule_id[iat] = imol
                        atom_type[iat] = int(temp[2]) - 1

                        for _ in range(1, num_atoms):
                            line = input_file.readline()
                            temp = line.split()
                            iat = int(temp[0]) - 1
                            imol = int(temp[1]) - 1
                            molecule_id[iat] = imol
                            atom_type[iat] = int(temp[2]) - 1

                        # rewind the LAMMPS data file
                        input_file.seek(0)
                        return

def read_bonds(input_file, num_bonds, sorted_bonds):
    """ store the atom-molecule and atom type info in two lists"""
    left_bonds, right_bonds = zeros(num_bonds, dtype=int), zeros(num_bonds, dtype=int)

    while True:
        line = input_file.readline()
        temp = line.split()
        if temp:
            if temp[0] == 'Bonds':
                while True:
                    line = input_file.readline()
                    temp = line.split()
                    if temp:
                        indx = int(temp[0]) - 1
                        left_bonds[indx] = int(temp[2]) - 1
                        right_bonds[indx] = int(temp[3]) - 1
                        for _ in range(1, num_bonds):
                            line = input_file.readline()
                            temp = line.split()
                            indx = int(temp[0]) - 1
                            left_bonds[indx] = int(temp[2]) - 1
                            right_bonds[indx] = int(temp[3]) - 1
                        num_atoms = max(amax(left_bonds), amax(right_bonds)) + 1
                        for iat in range(0, num_atoms):
                            sorted_bonds.append(concatenate((left_bonds[right_bonds == iat], \
                                             right_bonds[left_bonds == iat])))
                        input_file.seek(0) # rewind the LAMMPS data file
                        return
