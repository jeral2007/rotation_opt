import scipy as sc


def parse_line(line):
    if '$end' in line:
        return ()
    aux = line.split()
    aux[0:3] = [float(x) for x in aux[0:3]]
    if len(aux) < 5:
        return (aux, 'FROZEN')
    else:
        return (aux, 'UNFROZEN')


def parse_pattern(filename):
    f = open(filename)
    line = f.next()
    if '$coord' not in line:
        raise ValueError

    frozen_ats = {}
    at_blocks = {}  # atom blocks storage
    a_s = []  # reference to the block of atom with given number, it index in it and name of the atom
    frozen_ats['c'] = sc.array([0e0, 0e0, 0e0])  # center of frozen atoms subsystem is arbitrary and deliberatily assigned to 0e0 0e0 0e0
    frozen_ats['coords'] = []
    frozen_ats['new_coords'] = []
    for line in f:
        if '$end' in line:
            break
        atom, atype = parse_line(line)
        if atype == 'FROZEN':
            frozen_ats['coords'] += [sc.array(atom[0:3])]
            frozen_ats['new_coords'] +=[sc.array(atom[0:3])]

            a_s += [(frozen_ats, len(frozen_ats['coords'])-1, atom[3])]
        elif '_b' in atom[4]:
            tmp = atom[4].split('b')
            bnum = int(tmp[1])
            if bnum not in at_blocks.keys():
                at_blocks[bnum] = {}
                at_blocks[bnum]['coords'] = []
                at_blocks[bnum]['new_coords'] = []
            if 'c' in atom[4]:
                at_blocks[bnum]['c'] = sc.array(atom[0:3])

            at_blocks[bnum]['coords'] += [sc.array(atom[0:3])]
            at_blocks[bnum]['new_coords'] += [sc.array(atom[0:3])]
            a_s += [(at_blocks[bnum],
                     len(at_blocks[bnum]['coords'])-1, atom[3])]
        else:
            raise ValueError

    for block, bi, _ in a_s:
        block['coords'][bi] -= block['c']
    for block in at_blocks.itervalues():
        block['grads'] = [sc.zeros(3) for c in block['coords']]

    frozen_ats['grads'] = [sc.zeros(3) for c in frozen_ats['coords']]
    f.close()
    return a_s, at_blocks, frozen_ats


def parse_control(a_s, filename):
    f = open(filename)
    for line in f:
        if "$grad          cartesian gradients" in line:
            break
    en = 0
    while(True):
        line = f.next()
        if '$end' in line:
            f.close()
            return a_s, en
        en = float(line.split('=')[2].split()[0])
        for at, line in zip(a_s, f):
            if at[2] != line.split()[3]:
                print at[2], line.split()[3]
        for at, line in zip(a_s, f):
            aux = line.replace('D', 'e').replace('d', 'e').split()[0:3]
            block, ind, _ = at
            block['grads'][ind] = sc.array(map(float, aux))


def write_coord(a_s, filename='coord'):
    f = open(filename, 'w')
    f.write('$coord\n')
    for block, bl_ind, name in a_s:
        f.write("     {0[0]:.12f} {0[1]:.12f} {0[2]:.12f}  {1}\n".format(
            block['new_coords'][bl_ind], name))
    f.write('$end\n')
    f.close()
