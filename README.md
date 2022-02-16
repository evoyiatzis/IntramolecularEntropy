This is a python code to compute the intramolecular entropy of a molecule using the Schlitter quasi-harmonic approach.

It is a serial code which makes use of NumPy for most of the numerical calculations. The input files have the same format as those obtained from LAMMPS software.

The software has been tested using Python 3.6 with Numpy 1.11.3. It can be executed by issuing the following command at the prompt:
	python schlitter.py temperature output_entropy_file input_data_file type_of_dump_files list_of_input_dump_files --restart {positive integer number} --units {string}
	
A short description of the code as well as a complete list of the input, optional and output arguments can be obtained at the command line by executing
	python schlitter.py --help
	
The required input parameters are the temperature of the system, in Kelvin, a LAMMPS data file, which specifies the connectivity of the particles in the system, a string specifying the format of the provided dump files to be analyzed and finally the name of a single or several dump files containing the trajectory. The structure of the LAMMPS data file should follow the one dictated by the “atom_style full” command. A support for xyz files as well as for native LAMMPS dump files with scaled or unscaled coordinates is provided. The dump files are expected to be of ASCII type as well as to be compressed with the zlib library. 
The output of the code is stored in two ASCII text files. The name of the first file is “ReferenceConfiguration.txt”; it is located in the same directory as the code. In this file, the configuration with the smallest radius of gyration for all molecules is stored. In a second file, which has a user specified name, the computed entropy of each molecule is saved followed by two columns containing the non-zero eigenvalues of the covariance matrix together with their contribution to intramolecular entropy. All entropy values are given in k_B units.
Two optional arguments can be specified. The first one is the keyword --units followed by a string which specifies the employed unit system in the simulation. If this keyword is not specified, then the code assumes that the unit system is “REAL”. The second optional argument is the keyword --restart followed by a positive integer number. If specified, the code assumes that the calculation is restarted. 
---

The code was developed at Technical University of Darmstadt.  It is an open-source code, distributed freely under the terms of the GNU Public License (GPL).

The latest version of the code is available at https://github.com/evoyiatzis/IntramolecularEntropy

The maintainer of the code is Evangelos Voyiatzis, who can be emailed at evoyiatzis@gmail.com.  

The distribution includes the following files and directory:

README.md: this file

LICENSE: the GNU open-source license

Source: directory with the source code

