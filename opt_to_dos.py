from __future__ import print_function
import sys
import subprocess
import itertools
import re


"""
Extract optimized geometry from the end of an opt log file
Combine with header and pseudohydrogen info from the same opt gjf file
Write to a new file
"""


def get_line_number(log_file, search_str):
    """
    Search log_file for search_str and return the line number of the first occurrence.
    """

    line = subprocess.Popen(['grep', '-m 1', '-n', search_str, log_file], stdout=subprocess.PIPE).communicate()[0].decode(sys.stdout.encoding)
    if line is None:
        raise ValueError("Search term " + search_str + "not found in log file.")
    return int(re.search('\d+', line).group(0))


def is_float(value):
    """
    Return true if input string can be converted to a float.
    """

    try:
        float(value)
        return True
    except ValueError:
        return False


def is_atom_line(line):
    """
    Atom lines are formatted:
    Cu    1.23456    -1.23456    6.54321
    """

    if len(line.split()) != 4:  # 4 things in line: element, x, y, z
        return False
    elif len(line.split()[0]) > 2:  # Element symbol is 1 or 2 characters
        return False
    elif not is_float(line.split()[1]) or \
         not is_float(line.split()[1]) or \
         not is_float(line.split()[1]):
        return False
    else:
        return True


opt_filename = sys.argv[1]
opt_log = opt_filename + ".log"
opt_gjf = opt_filename + ".gjf"

if len(sys.argv) > 2:
    out_filename = sys.argv[2]
else:
    out_filename = opt_filename + "_optgeom.gjf"

# Read all lines from input file into a list
with open(opt_gjf, 'rb') as f:
    gjf_lines = f.readlines()

# Read end of log file and put atom types/coordinates into a list
geom_start = get_line_number(opt_log, "l9999.exe")
log_lines = ""
with open(opt_log, 'rb') as f:
    for line in itertools.islice(f, geom_start, None):
        log_lines += line.strip()
log_atoms = log_lines.split('\\')


with open(out_filename, 'wb') as f:
    for line in gjf_lines:
        line = line.strip()
	if is_atom_line(line):
            gjf_atom_type = line.split()[0].strip()  # get the atom type

            while True:
                if log_atoms[0].split(',')[0] == gjf_atom_type:  # if output/input atom types match
                    break
                else:
                    log_atoms.pop(0) # remove and discard line from log file

            log_atom_list = log_atoms.pop(0).split(',')  # splits into atom type, x, y, z
            x = float(log_atom_list[1])
            y = float(log_atom_list[2])
            z = float(log_atom_list[3])

            f.write('{0:8}{1:12.6f}{2:12.6f}{3:12.6f} \n'.format(gjf_atom_type, x, y, z))
	
	elif line.startswith("%chk="):
	    chk_filename = out_filename.split('.')[0] + '.chk'
	    f.write('%chk=' + chk_filename + '\n')

        elif line.startswith("#p"):  # keyword line
	    old_keyword_list = line.split()
	    new_keyword_line = ""
	    for keyword in old_keyword_list:
		if "opt" in keyword:
		    new_keyword_line += "GFInput Iop(3/33=1) Pop(Full,NPA) "
		else:
		    new_keyword_line += (keyword + " ")
	    f.write(new_keyword_line.strip() + '\n')
	else:
            f.write(line + '\n')  # all other lines written from input file
