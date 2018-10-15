from __future__ import (absolute_import, division, print_function)

import itertools
import os.path
import re
import sys
import subprocess
import time
import numpy as np
import scipy.sparse

"""
getDOS.py runs in the command line on a gaussian log file.
Probably requires python 3.5 or 3.6.
Requires numpy and scipy.sparse.
"""


class Atom:
    """
    Stores information associated with an atom.

    atom_type : string, element identity (C, Cu, Se, etc.)
    atom_num : int, number of that atom
    AOs : list of integers corresponding to AO numbers associated with that atom
    """

    def __init__(self, atom_type, atom_num):
        self.atom_type = atom_type
        self.atom_num = atom_num
        self.AOs = []


class MolecularOrbital:
    """
    Stores information associated with a molecular orbital.

    coef_sums : dictionary
        keys : AO types (Cu_d, C_s, etc.)
        values : sum of contributions from that AO type
                 obtained from the product of the overlap and MO coef matrices
    energy : MO energy in Hartrees
    """

    def __init__(self, energy):
        self.coef_sums = {}
        self.energy = energy

    def setup_coef_sums(self, target_atom_type=None):
        """Initialize dictionary by reading keys from unique_atoms list"""
        if target_atom_type:  # only read one atom type
            for AO_type in unique_atoms[target_atom_type]:
                self.coef_sums[target_atom_type + "_" + AO_type] = 0
        else:  # read all atom types
            for atom in unique_atoms:
                for AO_type in unique_atoms[atom]:
                    self.coef_sums[atom + "_" + AO_type] = 0



def get_line_number(log_file, search_str):
    """
    Search log_file for search_str and return the line number of the first occurrence.
    """

    line = subprocess.check_output(['grep', '-m 1', '-n', search_str, log_file]).decode(sys.stdout.encoding)
    if line is None:
        raise ValueError("Search term not found in log file.")
    return int(re.search('\d+', line).group(0))


def get_line(log_file, search_str):
    """
    Search log_file for search_str and return the line containing the first occurrence
    """
    return subprocess.check_output(['grep', '-m 1', search_str, log_file]).decode(sys.stdout.encoding)


def get_mo_energy(mo_num):
    """
    Finds the energy of the input MO by searching the log file.
    log, start_coef, and nbasis must be defined already.

    :param mo_num: integer MO number
    :return: energy of that MO in eV
    """
    with open(log_filepath, 'r') as file:
        start_line = start_coef + (mo_num - 1) // 5 * (nbasis + 3) + 2
        for line in itertools.islice(file, start_line, start_line + 1):
            return float(line.split()[(mo_num - 1) % 5 + 2]) * 27.2114


def get_mo_num(mo_energy_h):
    """
    Finds the number of the MO closest to (just below) the input energy.
    log, start_coef, end_coef, and nbasis must be defined already.

    :param mo_energy_h: energy in Hartrees (float)
    :return: MO number (int)
    """
    with open(log_filepath, 'r') as file:
        start_line = start_coef
        end_line = end_coef
        current_mo_num = 0
        for line in itertools.islice(file, start_line + 2, end_line, nbasis + 3):
            line = re.findall(r'-?[0-9]+\.[0-9]+', line)
            for eigenval in line:
                current_mo_num += 1
                if float(eigenval) > mo_energy_h:
                    return current_mo_num
        return current_mo_num


def get_arg():
    """
    Pops and returns the third argument in sys.argv.
    sys.argv should be ['getDOS.py', '-flags', arguments of flags, txt or log, txt or log]
    so this function will access the arguments of the flags, in order.
    Exits if there is not a third argument.
    """
    if len(sys.argv) > 2:
        return sys.argv.pop(2)
    else:
        print("Missing arguments.")
        sys.exit()


def write_output(file_object, MO_list, extra_name=''):
    """
    Print formatted output from a list of MolecularOrbital objects.
    Print headers from coef_sums keys, and data from coef_sums values for each MO.
    """

    file_object.write('{}\t{:8}\t'.format('MO num', 'Energy'))
    for AO_type in sorted(MO_list[0].coef_sums.keys()):
        if extra_name: # insert extra text after atom type (Cu1_s, Cu1_p, etc.)
            temp = AO_type.split("_")
            AO_type = temp[0] + extra_name + "_" + temp[1]
        file_object.write('{:8}\t'.format(AO_type))
    file_object.write('{:8}\t'.format('Occupied'))
    file_object.write('\n')

    # Data
    for MO_idx, MO in enumerate(MO_list):
        file_object.write(str(MO_idx + start_mo_num) + '\t')
        file_object.write('{:8.4f} \t'.format(MO.energy))  # 8 characters, 4 after decimal
        for AO_type in sorted(MO.coef_sums.keys()):
            file_object.write('{:8.4f} \t'.format(MO.coef_sums[AO_type]))

        if MO_idx + start_mo_num <= num_electrons:  # occupied or unoccupied
            file_object.write('{:8}'.format('1'))
        else:
            file_object.write('{:8}'.format('0'))

        file_object.write('\n')






def print_help():
    """Print help"""
    print()
    print("Example: python getDOS.py -e -10 10 logfilename.log outputfilename.txt")
    print("\t-e for energy range")
    print("\t-n for MO numbers")
    print("\t-i to get MO information only")
    print("\t-b to get beta MOs for open-shell systems")
    print("\t-a to output contributions from a specific atom's AOs")
    print("With no flags, default is alpha HOMO +/- 10 MOs.")
    print("If no output filename is specified, logfilename_dos.txt is used.")


################################################
### Processing user input and checking files ###
################################################
start_time = time.time()

# Example sys.argv ['getDOS.py', '-e', '-10', '10', 'filename.log', 'filename.txt']
# or [script name, flags, flag arguments, log or text filenames]

if len(sys.argv) < 2:
    print_help()
    sys.exit("Please enter a log file.")

# Get flags if the user entered any
flag = ''
if len(sys.argv) > 2:
    flag = sys.argv[1]

# Get log filename and optional output filename from user
name1, ext1 = os.path.splitext(sys.argv[-1])  # check the last input argument
if ext1 == ".log":
    log_filepath = sys.argv.pop()
elif ext1 == ".txt":
    out_filepath = sys.argv.pop()
else:
    sys.exit("Last argument is neither a log file nor a text file.")

name2, ext2 = os.path.splitext(sys.argv[-1])  # check the second-to-last (now last) input argument
if ext2 == ".log":
    log_filepath = sys.argv.pop()
elif ext2 == ".txt":
    out_filepath = sys.argv.pop()

try:  # if log is not assigned yet, error
    log_filepath
except NameError:
    sys.exit("Couldn't find log file in input.")

try:  # if outfile is not assigned yet, set it to logfilename_getDOS.txt
    out_filepath
except NameError:
    log_name, log_ext = os.path.splitext(log_filepath)
    out_filepath = log_name + "_getDOS.txt"

# Make sure log file exists
if not os.path.isfile(log_filepath):
    print("Log file {} does not exist.".format(log_filepath))
    sys.exit()

# -------------------------#

# Get number of alpha electrons from log file
line = get_line(log_filepath, 'alpha electrons')  # Line: "5 alpha electrons 5 beta electrons"
(num_alp, num_bet) = re.findall(r'\d+', line)  # get numbers in this line
num_alp = int(num_alp)
num_bet = int(num_bet)

if num_alp == num_bet: # closed-shell system
    beta = False
    mo_type = ''
    num_electrons = num_alp
else: # open-shell system
    if 'b' in flag: # user specified beta
        beta = True
        mo_type = 'beta '
        num_electrons = num_bet
        print("Open-shell system, finding beta MOs")
    else: # no specification, so do alpha
        beta = False
        mo_type = 'alpha '
        num_electrons = num_alp
        print("Open-shell system, finding alpha MOs")

# Get number of basis functions from log file
line = get_line(log_filepath, 'NBasis')
nbasis = int(re.search('\d+', line).group(0)) # Get the first number from this line

# Find the lines of the log file containing the overlap matrix
start_overlap = get_line_number(log_filepath, '*** Overlap ***') + 1
x = -(-nbasis // 5)
end_overlap = start_overlap + x + nbasis * x - ((x-1)*x//2*5) - 1
#print("Overlap matrix is from line {} to {}".format(start_overlap, end_overlap))

# Find the lines of the log file containing the MO coefficient matrix
if beta:
    start_coef = get_line_number(log_filepath, 'Beta Molecular Orbital Coefficients')
else:
    start_coef = get_line_number(log_filepath, 'Orbital Coefficients')
end_coef = start_coef + (-(-nbasis // 5) * (nbasis + 3))
# print("MO coef matrix is from line {} to {}".format(start_coef, end_coef))

#-------------------------------------------#

# Handle other flags in user input (i, e, n, a)

# Print info and exit if the user specified -i (information only)
if 'i' in flag:
    print("Number of alpha electrons: {}".format(num_alp))
    print("Number of beta electrons: {}".format(num_bet))
    print("Number of basis functions: {}".format(nbasis))
    print("There are {0} {1}MOs and the {1}HOMO is MO {2}.".format(nbasis, mo_type, num_electrons))
    print("The {}HOMO energy is {:6.2f} eV.".format(mo_type, get_mo_energy(num_electrons)))
    print("The {}LUMO energy is {:6.2f} eV.".format(mo_type, get_mo_energy(num_electrons + 1)))
    sys.exit()


# Handle other flags in order
do_atom = False
for letter in flag:
    if letter == 'e':  # specify energy range
        start_en = float(get_arg())
        start_mo_num = get_mo_num(start_en / 27.2114)
        end_en = float(get_arg())
        end_mo_num = get_mo_num(end_en / 27.2114)

    elif letter == 'n':  # specify number of atoms
        start_mo_num = int(get_arg())
        end_mo_num = int(get_arg())
        if start_mo_num < 1 or start_mo_num > nbasis:
            sys.exit("{} is not a valid starting MO number.".format(start_mo_num))
        if end_mo_num > nbasis:
            sys.exit("{} is not a valid ending MO number. \n"
                     "Ending MO number should be {} or smaller.".format(end_mo_num, nbasis))
        if start_mo_num > end_mo_num:
            sys.exit("Starting MO number should not be greater than ending MO number.")

    elif letter == 'a': # set target atom
        do_atom = True
        target_atom_num = int(get_arg())
        if target_atom_num < 1:
            sys.exit("{} is not a valid atom number.".format(target_atom_num))
        print("Calculating individual contributions for atom {}.".format(target_atom_num))

# if MO numbers are not defined yet, set them to HOMO +/- 10
try:
    start_mo_num
except:
    start_mo_num = max(num_electrons - 10, 1)
try:
    end_mo_num
except:
    end_mo_num = min(num_electrons + 10, nbasis)


print("Reading data from {}".format(os.path.basename(log_filepath)))
print("Writing output to {}".format(os.path.basename(out_filepath)))
print("There are {0} {1}MOs and the {1}HOMO is MO {2}.".format(nbasis, mo_type, num_electrons))
print("Calculating {}MOs {} to {}".format(mo_type, start_mo_num, end_mo_num))
print()


# -------------------------#

# Build atoms, unique_atoms, and AOs from MO coef matrix section of log file
atoms = []  # list of atom objects
unique_atoms = {}  # dictionary mapping 'Cu' --> ['s', 'p', 'd'] etc.
AOs = []  # list of all AO types Cu_s, Cu_p, etc. with index corresponding to AO num - 1

with open(log_filepath, 'r') as file:
    for line in itertools.islice(file, start_coef + 3, start_coef + nbasis + 3): # use first piece of MO coef matrix
        line = line[:20].split()
        ao_num = int(line[0])

        if len(line) == 4:  # example:  1 1  C  1S (AO number, atom number, atom type, AO type)
            add_unique_aos = False  # default: don't add AOs to the unique atom list
            atom_num = int(line[1])
            atom_type = line[2]
            ao_type = line[3][1].lower()
            atoms.append(Atom(atom_type, atom_num))  # add this atom to the list of Atom objects

            if atom_type not in unique_atoms:  # add to list of unique atoms
                unique_atoms[atom_type] = []
                unique_atoms[atom_type].append(ao_type)
                add_unique_aos = True  # continue adding AOs for this atom type

        else:  # line contains AOs only (AO number, AO type)
            if add_unique_aos and ao_type not in unique_atoms[atom_type]:
                unique_atoms[atom_type].append(ao_type)
            ao_type = line[1][1].lower()

        atoms[atom_num - 1].AOs.append(ao_num)  # add this AO number to the atom's AO list
        AOs.append(atom_type + "_" + ao_type)  # add this AO identity (Cu_s) to the list of AOs

# -------------------------#

# Read the overlap matrix from the log file
start_overlap_time = time.time()
with open(log_filepath, 'r') as file:
    overlap_matr = scipy.sparse.lil_matrix((nbasis, nbasis))
    for line in itertools.islice(file, start_overlap - 1, end_overlap): # -1 to get index from line num
        if 'D' not in line:  # line contains integers (section header: basis function numbers)
            col_j0 = int(line.split()[0]) - 1  # first integer gives column index
        else:  # line contains numbers in scientific notation (overlap values)
            col_j = col_j0  # reset column counter
            line = line.replace('D', 'E').split()  # change scientific notation from D to E
            row_i = int(line[0]) - 1  # first integer gives row index
            for element in line[1:]:
                overlap_matr[row_i, col_j] = element
                if row_i != col_j:  # off diagonal
                    overlap_matr[col_j, row_i] = element # set symmetric element too
                col_j += 1  # increment column

end_overlap_time = time.time()
print("Read overlap into matrix in {:6.3f} seconds".format(end_overlap_time - start_overlap_time))

# -------------------------#

# Read the MO coefficient matrix and eigenvalues from the log file
start_mocoef_time = time.time()
with open(log_filepath, 'r') as file:
    mocoef_matr = np.empty((nbasis, end_mo_num - start_mo_num + 1))
    eigenvalues = []

    # MO coef matrix part of log file is a repeating pattern of lines
    # so we can find the exact lines that contain the specified MOs
    start_line = start_coef + 1 + (start_mo_num - 1) // 5 * (nbasis + 3)
    end_line = start_coef + 1 + ((end_mo_num - 1) // 5 + 1) * (nbasis + 3)

    for line in itertools.islice(file, start_line - 1, end_line - 1):
        if line[0:21].isspace():  # first 21 chars are blank
            if line.split()[0].isdigit():  # line contains numbers (MO numbers)
                col_j0 = int(line.split()[0]) - start_mo_num  # use to get column index
        elif "Eigenvalues" in line:  # line containing eigenvalues
            col_j = col_j0  # reset column counter
            for eigenval in line[21:].split():
                if 0 <= col_j <= end_mo_num - start_mo_num: # if column is in desired MO range
                    eigenvalues.append(float(eigenval))
                col_j += 1 # increment column
        else:  # line containing MO coef values
            col_j = col_j0  # reset column counter
            row_i = int(line.split()[0]) - 1  # first integer gives row index
            for idx, element in enumerate(line[21:].split()):
                if 0 <= col_j <= end_mo_num - start_mo_num: # if column is in desired MO range
                    mocoef_matr[row_i, col_j] = element
                col_j += 1  # increment column
end_mocoef_time = time.time()
print("Read MO coefs into matrix in {:6.3f} seconds".format(end_mocoef_time - start_mocoef_time))

# Set up list of MOs from eigenvalues
MOs = []  # list of MO objects
for eigenval in eigenvalues:
    new_MO = MolecularOrbital(float(eigenval) * 27.2114)
    new_MO.setup_coef_sums() # initialize dictionary
    MOs.append(new_MO)

# -------------------------#

start_matrix_time = time.time()
# Do matrix multiplication
# coefs_T * overlap, product multiplied elementwise by coefs_T
result = np.multiply((mocoef_matr.transpose() * overlap_matr.tocsc()), (mocoef_matr.transpose()))
end_matrix_time = time.time()
print("Did matrix multiplication in {:6.3f} seconds".format(end_matrix_time - start_matrix_time))

# -------------------------#

# Current data structures
# -----------------------
# atoms : list of Atom objects
#       each Atom contains .atom_num, .atom_type, .AOs (list of AO numbers)
# unique_atoms : dictionary mapping atom type to list of AO types
#       example: Cu --> [s, p, d]
# MOs : list of MolecularOrbital objects
#       each MolecularOrbital contains .energy and .coef_sums
#       coef_sums is a dictionary mapping each atom and AO type to a value (currently 0)
#           example: Cu_d --> 0.000
# result : a matrix containing AO contributions to MOs including overlap
#       each row is a MO
#       element i, j is the contribution of AO j to MO (i + start_mo_num)

# -------------------------#

# Read data from result matrix into MO coef_sums
start_sum_time = time.time()
for MO_idx, MO in enumerate(MOs):
    result_row = result[MO_idx]  # each row of the result matrix is one MO
    for AO_idx, AO_type in enumerate(AOs):
        MO.coef_sums[AO_type] += 100 * result_row[AO_idx]  # each column of the result matrix is one AO
end_sum_time = time.time()
print("Summed overlap contributions in {:6.3f} seconds".format(end_sum_time - start_sum_time))

# -------------------------#

# Format and print output
with open(out_filepath, 'w') as f:
    # Headers
    f.write('{}\t{:8}\t'.format('MO num', 'Energy'))
    for AO_type in sorted(MOs[0].coef_sums.keys()):
        f.write('{:8}\t'.format(AO_type))
    f.write('{:8}\t'.format('Occupied'))
    f.write('\n')

    # Data
    for MO_idx, MO in enumerate(MOs):
        f.write(str(MO_idx + start_mo_num) + '\t')
        f.write('{:8.4f} \t'.format(MO.energy))  # 8 characters, 4 after decimal
        for AO_type in sorted(MO.coef_sums.keys()):
            f.write('{:8.4f} \t'.format(MO.coef_sums[AO_type]))

        if MO_idx + start_mo_num <= num_electrons:  # occupied or unoccupied
            f.write('{:8}'.format('1'))
        else:
            f.write('{:8}'.format('0'))

        f.write('\n')

print("Done!")

# -------------------------#

# Calculate individual atom contributions if needed
print()
while do_atom:
    if target_atom_num - 1 > len(atoms):
        sys.exit("Requested atom number is out of range.")
    target_atom_type = atoms[target_atom_num - 1].atom_type
    print("Calculating contribution from atom {} {}.".format(target_atom_type, target_atom_num))

    start_AO = atoms[target_atom_num - 1].AOs[0] - 1 # index
    end_AO = atoms[target_atom_num - 1].AOs[-1] # index, exclusive
    #print("Using AOs {} to {}".format(start_AO + 1, end_AO))

    target_atom_MOs = []  # list of MO objects but ONLY storing contributions from the target atom
    for eigenval in eigenvalues:
        new_MO = MolecularOrbital(float(eigenval) * 27.2114)
        new_MO.setup_coef_sums(target_atom_type)
        target_atom_MOs.append(new_MO)

    for MO_idx, MO in enumerate(target_atom_MOs):
        result_row = result[MO_idx]  # each row of the result matrix is one MO
        for AO_idx in range(start_AO, end_AO):
            AO_type = AOs[AO_idx]
            MO.coef_sums[AO_type] += 100 * result_row[AO_idx]  # each column of the result matrix is one AO

    out_atom_filepath = os.path.splitext(out_filepath)[0] + "_" + target_atom_type + str(target_atom_num) +".txt"
    with open(out_atom_filepath, 'w') as f:
        write_output(f, target_atom_MOs, str(target_atom_num))

    if (time.time() - start_time) < 1: # if the calculation was fast, done
        break
    else: # if the calculation took some time, option to do more atoms
        next_atom = input("Do another atom? Enter atom number, return to do next atom, or any other input to stop: ")
        if next_atom == '':
            target_atom_num = target_atom_num + 1
        elif not next_atom.isnumeric():
            break
        else:
            target_atom_num = int(next_atom)


#--------------------------#

end_time = time.time()
print()
print("Total time: {:6.2f} seconds".format(end_time - start_time))
