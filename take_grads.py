#!/usr/bin/python 
from math import sqrt
import sys

def get_components(line):
  fmt = '-0.1274010D-01'
  i = 0
  as1 = ''
  res = []
  for c in line[7:]:
     if c =='D' or c == 'd':
         as1 += 'e'
     else:
         as1 += c
     i+=1
     if i == len(fmt):
         res += [float(as1)]
         i = 0
         as1 = ''
  return res


def usage():
    usage_str = """Evaluates square root of the sum of atom  gradient squares, except of the atoms with names in atom_lst. 
Usage: {} grad_file coord_file [atom_lst]\n""".format(sys.argv[0])
    sys.stderr.write(usage_str)


def make_exclude_lst(coord_filename, lst):
    if len(lst) == 0 or lst is None:
       usage()
       quit()
    try:
        coord_file = open(coord_filename, 'r')
    except Exception:
        usage()
        quit()
    line = coord_file.next()
    if '$coord' not in line:
        sys.stderr.write('Incorrect coord file {}'.format(coord_filename))
        quit()
    result = []
    at_names = []
    for at_num, line in enumerate(coord_file):
        if '$end' in line:
            break
        aux = line.split()[3]
	if aux.lower() in lst:
            result +=[at_num + 1]
        else:
            at_names += [aux.lower()]
    return set(result), at_names


try:
    assert(len(sys.argv)>2)
    grad_file = open(sys.argv[1], 'r')
except Exception:
    usage()
    quit()
if len(sys.argv)>3:
    exclude_atoms = set(atom.lower() for atom in sys.argv[3:])
    exclude_lst, at_names = make_exclude_lst(sys.argv[2], exclude_atoms)
else:
    exclude_lst, at_names = make_exclude_lst(sys.argv[2], [])

for line in grad_file:
   if 'cartesian gradient of the energy' in line:
       break

grad_file.next()
grad_file.next()

print (exclude_lst)
print ('-'*36)
xg, yg, zg = [], [], []
while(True):
   line = grad_file.next()
   if 'resulting FORCE' in line:
       break
   aux = line.split()
   numbers = map(int, aux[1::2])
   line = grad_file.next()
   aux = get_components(line)
   xg += [float(t) for t, n in zip(aux, numbers) 
	     if n not in exclude_lst]
   line = grad_file.next()
   aux = get_components(line)
   yg += [float(t) for t, n in zip(aux, numbers) 
	     if n not in exclude_lst]
   line = grad_file.next()
   aux = get_components(line)
   zg += [float(t) for t, n in zip(aux, numbers) 
	     if n not in exclude_lst]
   line = grad_file.next()

for x, y, z, at in zip(xg, yg, zg, at_names):
    print "{0:>+2.7f} {1:>+2.7f} {2:>+2.7f} {3: >3}".format(x, y, z, at)
print "-"*36
print(sum(x**2+y**2+z**2 for x,y,z in zip(xg,yg,zg)))
print(sqrt(sum(x**2+y**2+z**2 for x,y,z in zip(xg,yg,zg))))
