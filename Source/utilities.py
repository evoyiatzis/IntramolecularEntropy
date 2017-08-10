#!/usr/local/bin/python

""" This module contains the functions:
- set_units
- kabsch
- form_displacement_matrix
"""

import sys
import numpy

def set_units(unit_system):
    """This function takes care of the units in the boltzmann and planck constant """

    constants = {"SI": ((6.62606957 * pow(10, -34))/(2*numpy.pi), 1.3806488 * pow(10, -23)), \
    "LJ": ((0.18292026)/(2*numpy.pi), 1.0), \
    "METAL": ((1.3807 * pow(10, -16))/(2*numpy.pi), 8.6173 * pow(10, -5)), \
    "CGS": ((6.6261 * pow(10, -27))/(2*numpy.pi), 1.3807 * pow(10, -16)), \
    "REAL": ((39.9021909505 * pow(10, -3))/(2*numpy.pi), 8.3142670736 * pow(10, -7)), \
    "NANO": (6.62606957 * pow(10, -4)/(2*numpy.pi), 1.3806488 * pow(10, -2)), \
    "MICRO": ((6.62606957 * pow(10, -13))/(2*numpy.pi), 1.3806488 * pow(10, -8)), \
    "ELECTRON": (0.142534616294 /(2*numpy.pi), 2.96993934136 * pow(10, -6))}

    if unit_system in constants:
        return constants[unit_system]
    else:
        sys.exit('This unit system does not exist !')

def kabsch(atomic_coord, ref):
    """It applies the kabsch algorithm to determine the rotated coordinates of a molecule """

    p_matrix = numpy.column_stack((atomic_coord[:, 0], atomic_coord[:, 1], atomic_coord[:, 2]))
    q_matrix = numpy.column_stack((ref[:, 0], ref[:, 1], ref[:, 2]))
    cov_matrix = numpy.dot(numpy.transpose(p_matrix), q_matrix)
    v_svd, s_svd, w_svd = numpy.linalg.svd(cov_matrix)
    if (numpy.linalg.det(v_svd) * numpy.linalg.det(w_svd)) < 0.0:
        s_svd[-1] = -s_svd[-1]
        v_svd[:, -1] = -v_svd[:, -1]

    xopt2 = numpy.dot(v_svd, w_svd)
    p_all = numpy.dot(p_matrix, xopt2)

    return p_all

def form_disp_matrix(atomic_coord, square_mass):
    """"This function computes the displacement matrix of a single configuration"""
    atom_coords = numpy.concatenate([atomic_coord[:, 0], atomic_coord[:, 1], atomic_coord[:, 2]])
    sq_ext_masses = numpy.concatenate([square_mass, square_mass, square_mass])
    disp_matrix = numpy.outer(sq_ext_masses, sq_ext_masses)*numpy.outer(atom_coords, atom_coords)

    return disp_matrix
