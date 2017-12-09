#!/usr/local/bin/python

"""This file contains the functions
to read from the dump files
"""

from itertools import islice

def dump_file(input_file, num_atoms, atom_coord):
    """This function reads configurations from a LAMMPS dump file
     where coordinates are stored in unscaled format"""

    lines_gen = islice(input_file, num_atoms+9)

    icount = 0
    for line in lines_gen:
        icount += 1
        if icount == 6:
            temp = line.split()
            xlength = float(temp[1]) - float(temp[0])
        if icount == 7:
            temp = line.split()
            ylength = float(temp[1]) - float(temp[0])
        if icount == 8:
            temp = line.split()
            zlength = float(temp[1]) - float(temp[0])
        if icount < 10:
            continue
        temp = line.split()
        indx = int(temp[0]) - 1
        if len(temp) == 7:
            atom_coord[indx, 0] = float(temp[1]) + int(temp[4]) * xlength
            atom_coord[indx, 1] = float(temp[2]) + int(temp[5]) * ylength
            atom_coord[indx, 2] = float(temp[3]) + int(temp[6]) * zlength
        else:
            atom_coord[indx, 0] = float(temp[1])
            atom_coord[indx, 1] = float(temp[2])
            atom_coord[indx, 2] = float(temp[3])

    if icount == (num_atoms+9):
        return False

    return True

def scaled_dump_file(input_file, num_atoms, atom_coord):
    """This function reads configurations from a LAMMPS dump file where
     coordinates are stored in scaled format"""

    lines_gen = islice(input_file, num_atoms+9)

    icount = 0
    for line in lines_gen:
        icount += 1
        if icount == 6:
            temp = line.split()
            xlength = float(temp[1]) - float(temp[0])
        if icount == 7:
            temp = line.split()
            ylength = float(temp[1]) - float(temp[0])
        if icount == 8:
            temp = line.split()
            zlength = float(temp[1]) - float(temp[0])
        if icount < 10:
            continue
        temp = line.split()
        indx = int(temp[0]) - 1
        if len(temp) == 7:
            atom_coord[indx, 0] = float(temp[1]) * xlength + int(temp[4]) * xlength
            atom_coord[indx, 1] = float(temp[2]) * ylength + int(temp[5]) * ylength
            atom_coord[indx, 2] = float(temp[3]) * zlength + int(temp[6]) * zlength
        else:
            atom_coord[indx, 0] = float(temp[1]) * xlength
            atom_coord[indx, 1] = float(temp[2]) * ylength
            atom_coord[indx, 2] = float(temp[3]) * zlength

    if icount == (num_atoms+9):
        return False

    return True

def xyz_file(input_file, num_atoms, atom_coord):
    """
This function reads input configuration with a xyz file format.
    """
    lines_gen = islice(input_file, num_atoms+9)

    iat, icount = 0, 0
    for line in lines_gen:
        icount += 1
        if icount < 2:
            continue
        temp = line.split()
        if len(temp) == 0:
            return True
        atom_coord[iat, 0] = float(temp[1])
        atom_coord[iat, 1] = float(temp[2])
        atom_coord[iat, 2] = float(temp[3])
        iat += 1

    if icount == (num_atoms+9):
        return False

    return True

def intermediate_file(input_file, num_atoms, atom_coord):
    """This function reads data from the intermediate file created by the code."""

    lines_gen = islice(input_file, num_atoms)

    iat = 0
    for line in lines_gen:
        temp = line.split()
        if len(temp) == 0:
            return True
        atom_coord[iat, 0] = float(temp[0])
        atom_coord[iat, 1] = float(temp[1])
        atom_coord[iat, 2] = float(temp[2])
        iat += 1

    if iat == num_atoms:
        return False

    return True

def read_configuration(string, input_file, num_atoms, atom_coord):
    """This function manages reading from files. It provides an interface
    to the other functions contained in the same module."""

    dict_func = {
        'xyz': xyz_file,
        'scaled': scaled_dump_file,
        'typical': dump_file,
        'intermediate': intermediate_file,
    }

    if string in dict_func:
        return dict_func[string](input_file, num_atoms, atom_coord)
    else:
        sys.exit('This type of dump file does not exist !')

