'''
Crystallographic symmetry related commands.

(c) 2010-2012 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def cellbasis(angles, edges):
    '''
DESCRIPTION

    API only. For the unit cell with given angles and edge lengths calculate
    the basis transformation (vectors) as a 4x4 numpy.array
    '''
    from math import cos, sin, radians, sqrt
    import numpy

    rad = [radians(i) for i in angles]
    basis = numpy.identity(4)
    basis[0][1] = cos(rad[2])
    basis[1][1] = sin(rad[2])
    basis[0][2] = cos(rad[1])
    basis[1][2] = (cos(rad[0]) - basis[0][1]*basis[0][2])/basis[1][1]
    basis[2][2] = sqrt(1 - basis[0][2]**2 - basis[1][2]**2)
    edges.append(1.0)
    return basis * edges # numpy.array multiplication!

def supercell(a=1, b=1, c=1, object=None, color='green', name='supercell',
        withmates=1, prefix='m'):
    '''
DESCRIPTION

    Draw a supercell, as requested by Nicolas Bock on the pymol-users
    mailing list (Subject: [PyMOL] feature request: supercell construction
    Date: 04/12/2010 10:12:17 PM (Mon, 12 Apr 2010 14:12:17 -0600))

USAGE

    supercell a, b, c [, object [, color [, name [, withmates]]]]

ARGUMENTS

    a, b, c = integer: repeat cell in x,y,z direction a,b,c times
    {default: 1,1,1}

    object = string: name of object to take cell definition from

    color = string: color of cell {default: blue}

    name = string: name of the cgo object to create {default: supercell}

    withmates = bool: also create symmetry mates in displayed cells
    {default: 1}

SEE ALSO

    show cell

    '''
    import numpy
    from pymol import cgo

    if object is None:
        object = cmd.get_object_list()[0]
    withmates = int(withmates)

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]

    basis = cellbasis(cell_angles, cell_edges)
    assert isinstance(basis, numpy.ndarray)

    ts = list()
    for i in range(int(a)):
        for j in range(int(b)):
            for k in range(int(c)):
                ts.append([i,j,k])

    obj = [
        cgo.BEGIN,
        cgo.LINES,
    ]

    groupname_fmt = prefix + ('%d_%d_%d' if any(int(i) > 9
            for i in [a, b, c]) else '%d%d%d')

    for t in ts:
        shift = basis[0:3,0:3] * t
        shift = shift[:,0] + shift[:,1] + shift[:,2]

        for i in range(3):
            vi = basis[0:3,i]
            vj = [
                numpy.array([0.,0.,0.]),
                basis[0:3,(i+1)%3],
                basis[0:3,(i+2)%3],
                basis[0:3,(i+1)%3] + basis[0:3,(i+2)%3]
            ]
            for j in range(4):
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j]).tolist())
                obj.append(cgo.VERTEX)
                obj.extend((shift + vj[j] + vi).tolist())

        if withmates:
            groupname = groupname_fmt % tuple(t)
            symexpcell(groupname + '_', object, *t)
            cmd.group(groupname, groupname + '_*')

    obj.append(cgo.END)

    if name != 'none':
        cmd.delete(name)
        cmd.load_cgo(obj, name)
        cmd.color(color, name)

def symexpcell(prefix='mate', object=None, a=0, b=0, c=0):
    '''
DESCRIPTION

    Creates all symmetry-related objects for the specified object that
    occur with their bounding box center within the unit cell.

USAGE

    symexpcell prefix, object, [a, b, c]

ARGUMENTS

    prefix = string: prefix of new objects

    object = string: object for which to create symmetry mates

    a, b, c = integer: create neighboring cell {default: 0,0,0}

SEE ALSO

    symexp
    '''
    import numpy
    from pymol import xray

    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]
    spacegroup = sym[6]

    basis = cellbasis(cell_angles, cell_edges)
    basis = numpy.matrix(basis)

    extent = cmd.get_extent(object)
    center = sum(numpy.array(extent)) * 0.5
    center = numpy.matrix(center.tolist() + [1.0]).T
    center_cell = basis.I * center

    spacegroup = xray.space_group_map.get(spacegroup, spacegroup)

    i = 0
    matrices = xray.sg_sym_to_mat_list(spacegroup)
    for mat in matrices:
        i += 1

        mat = numpy.matrix(mat)
        shift = -numpy.floor(numpy.array(mat * center_cell)[0:3, 0])
        shift = shift.flatten().astype(int)
        shift += [a, b, c]
        mat[0:3, 3] += shift.reshape((3, 1))

        mat = basis * mat * basis.I
        mat_list = list(mat.flat)

        name = '%s%d' % (prefix, i)
        cmd.create(name, object)
        cmd.transform_object(name, mat_list, 0)

        cmd.set_title(name, 1, "{}_{}{}{}".format(i, *(shift + 5).tolist()))

        if len(matrices) > 1:
            cmd.color(i + 1, name)

def pdbremarks(filename):
    '''
DESCRIPTION

    API only. Read REMARK lines from PDB file. Return dictionary with
    remarkNum as key and list of lines as value.
    '''
    remarks = dict()
    if not cmd.is_string(filename):
        f = filename
    elif filename[-3:] == '.gz':
        import gzip
        f = gzip.open(filename)
    else:
        f = open(filename)
    for line in f:
        recname = line[0:6]
        if recname == 'REMARK':
            num = int(line[7:10])
            lstring = line[11:]
            remarks.setdefault(num, []).append(lstring)
    return remarks

def biomolecule(name=None, filename=None, prefix=None, number=1, suffix=None,
        quiet=0):
    '''
DESCRIPTION

    Create biological unit (quaternary structure) as annotated by the REMARK
    350 BIOMOLECULE record.

USAGE

    biomolecule name [, filename [, prefix [, number ]]]

ARGUMENTS

    name = string: name of object and basename of PDB file, if
    filename is not given {default: first loaded object}

    filename = string: file to read remarks from {default: <name>.pdb}

    prefix = string: prefix for new objects {default: <name>}

EXAMPLE

    fetch 1rmv, async=0
    biomolecule 1rmv
    '''
    import os
    try:
        from .importing import local_mirror_pdb
    except ImportError:
        local_mirror_pdb = lambda n: ""

    try:
        import numpy
    except ImportError:
        numpy = None

    number, quiet = int(number), int(quiet)

    if name is None:
        name = cmd.get_object_list()[0]
    if prefix is None:
        prefix = name
    if suffix is None:
        suffix = str(number)
    if filename is None:
        candidates = [
            '%s.pdb' % (name),
            '%s/%s.pdb' % (cmd.get('fetch_path'), name),
            local_mirror_pdb(name),
        ]
        for filename in candidates:
            if os.path.exists(filename):
                break
        else:
            raise CmdException('please provide filename')
        if not quiet:
            print('loading from %s' % (filename))

    remarks = pdbremarks(filename)
    if 350 not in remarks:
        raise CmdException('There is no REMARK 350 in ' + filename)

    current = 1
    biomt = {current: {}}
    chains = tuple()

    for line in remarks[350]:
        if line.startswith('BIOMOLECULE:'):
            current = int(line[12:])
            biomt[current] = {}
        elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
            chains = tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('                   AND CHAINS:'):
            chains += tuple(chain.strip() for chain in line[30:].split(','))
        elif line.startswith('  BIOMT'):
            row = int(line[7])
            num = int(line[8:12])
            vec = line[12:].split()
            vec = list(map(float, vec))
            biomt[current].setdefault(chains, dict()).setdefault(num, []).extend(vec)

    if number not in biomt or len(biomt[number]) == 0:
        raise CmdException('no BIOMOLECULE number %d' % (number))

    if numpy is not None:
        mat_source = numpy.reshape(cmd.get_object_matrix(name), (4,4))
        mat_source = numpy.matrix(mat_source)

    for chains, matrices in biomt[number].items():
        for num in matrices:
            mat = matrices[num][0:12]
            mat.extend([0,0,0,1])
            copy = '%s_%s_%d' % (prefix, suffix, num)
            if not quiet:
                print('creating %s' % (copy))
            cmd.create(copy, 'model %s and chain %s' % (name, '+'.join(chains)))
            cmd.alter(copy, 'segi="%d"' % (num))

            if numpy is not None:
                mat = mat_source * numpy.reshape(mat, (4,4)) * mat_source.I
                mat = list(mat.flat)

            cmd.transform_object(copy, mat)

    cmd.disable(name)
    cmd.group('%s_%s' % (prefix, suffix), '%s_%s_*' % (prefix, suffix))

from chempy import cpv

class PutCenterCallback(object):
    prev_v = None

    def __init__(self, name, corner=0):
        self.name = name
        self.corner = corner
        self.cb_name = cmd.get_unused_name('_cb')

    def load(self):
        cmd.load_callback(self, self.cb_name)

    def __call__(self):
        if self.name not in cmd.get_names('objects'):
            import threading
            threading.Thread(None, cmd.delete, args=(self.cb_name,)).start()
            return

        v = cmd.get_view()
        if v == self.prev_v:
            return
        self.prev_v = v

        t = v[12:15]

        if self.corner:
            vp = cmd.get_viewport()
            R_mc = [v[0:3], v[3:6], v[6:9]]
            off_c = [0.15 * v[11] * vp[0] / vp[1], 0.15 * v[11], 0.0]
            if self.corner in [2,3]:
                off_c[0] *= -1
            if self.corner in [3,4]:
                off_c[1] *= -1
            off_m = cpv.transform(R_mc, off_c)
            t = cpv.add(t, off_m)

        z = -v[11] / 30.0
        m = [z, 0, 0, 0, 0, z, 0, 0, 0, 0, z, 0, t[0] / z, t[1] / z, t[2] / z, 1]
        cmd.set_object_ttt(self.name, m)

def cell_axes(object=None, a_color='0xd95f02', b_color='0x1b9e77', c_color='0x7570b3', name='cell_axes'):
    '''
DESCRIPTION

    Draw arrows corresponding to the crystallographic axes.

USAGE

    cell_axes [ object, [, a_color [, b_color [, c_color, [ name]]]]

ARGUMENTS

    object = string: name of object to take cell definition from

    a_color = string: color of a-axis {default: '0xd95f02'}

    b_color = string: color of b-axis {default: '0x1b9e77'}

    c_color = string: color of c-axis {default: '0x7570b3'}

    name = string: name of the cgo object to create {default: cell_axes}

    '''
    from pymol import cgo

    if object is None:
        object = cmd.get_object_list()[0]

    sym = cmd.get_symmetry(object)
    cell_edges = sym[0:3]
    cell_angles = sym[3:6]

    basis = cellbasis(cell_angles, cell_edges)
    print(basis)


    cmd.set('auto_zoom', 0)

    w = 0.06 # cylinder width
    l = 0.75 # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter
    obj = []
    
    import numpy as np

    A,B,C = basis[:3,:3].T
    r = A / np.linalg.norm(A) 
    rgb = cmd.get_color_tuple(a_color)
    obj.extend([
        cgo.CYLINDER, 0.0, 0.0, 0.0, l*r[0], l*r[1], l*r[2], w, *rgb, *rgb,
        cgo.CONE, l*r[0], l*r[1], l*r[2], (h+l)*r[0], (h+l)*r[1], (h+l)*r[2], d, 0.0, *rgb, *rgb, 1.0, 1.0
    ])

    r = B / np.linalg.norm(B) 
    rgb = cmd.get_color_tuple(b_color)
    obj.extend([
        cgo.CYLINDER, 0.0, 0.0, 0.0, l*r[0], l*r[1], l*r[2], w, *rgb, *rgb,
        cgo.CONE, l*r[0], l*r[1], l*r[2], (h+l)*r[0], (h+l)*r[1], (h+l)*r[2], d, 0.0, *rgb, *rgb, 1.0, 1.0
    ])

    r = C / np.linalg.norm(C) 
    rgb = cmd.get_color_tuple(c_color)
    obj.extend([
        cgo.CYLINDER, 0.0, 0.0, 0.0, l*r[0], l*r[1], l*r[2], w, *rgb, *rgb,
        cgo.CONE, l*r[0], l*r[1], l*r[2], (h+l)*r[0], (h+l)*r[1], (h+l)*r[2], d, 0.0, *rgb, *rgb, 1.0, 1.0
    ])

    PutCenterCallback(name, 1).load()
    cmd.load_cgo(obj, name)

cmd.extend('cell_axes', cell_axes)
cmd.extend('supercell', supercell)
cmd.extend('symexpcell', symexpcell)
cmd.extend('biomolecule', biomolecule)

# tab-completion of arguments
cmd.auto_arg[0]['biomolecule'] = cmd.auto_arg[0]['pseudoatom']
cmd.auto_arg[3]['supercell'] = cmd.auto_arg[0]['pseudoatom']

