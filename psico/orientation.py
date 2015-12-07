'''
Orientation, displacement and angle measurments of helices and domains.

(c) 2010-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException
from chempy import cpv

if cmd.get_version()[1] < 1.2:
    # legacy support to run this script standalone with PyMOL 0.99
    def get_unused_name(name):
        import random
        return name + '%04d' % random.randint(0, 1000)
    STATE = 1
else:
    from pymol.cmd import get_unused_name
    STATE = -1

def _vec_sum(vec_list):
    # this is the same as
    # return numpy.array(vec_list).sum(0).tolist()
    vec = cpv.get_null()
    for x in vec_list:
        vec = cpv.add(vec, x)
    return vec

def _common_orientation(selection, center, vec, visualize=1, scale=1.0, quiet=1):
    if visualize:
        visualize_orientation(vec, center, scale, True)
        cmd.zoom(selection, buffer=2)
    if not quiet:
        print(' Center: (%.2f, %.2f, %.2f) Direction: (%.2f, %.2f, %.2f)' % tuple(center + vec))

def visualize_orientation(direction, center=[0.0]*3, scale=1.0, symmetric=False, color='green', color2='red'):
    '''
DESCRIPTION

    Draw an arrow. Helper function for "helix_orientation" etc.
    '''
    from pymol import cgo

    color_list = cmd.get_color_tuple(color)
    color2_list = cmd.get_color_tuple(color2)

    if symmetric:
        scale *= 0.5
    end = cpv.add(center, cpv.scale(direction, scale))
    radius = 0.3

    obj = [cgo.SAUSAGE]
    obj.extend(center)
    obj.extend(end)
    obj.extend([
        radius,
        0.8, 0.8, 0.8,
    ])
    obj.extend(color_list)

    if symmetric:
        start = cpv.sub(center, cpv.scale(direction, scale))
        obj.append(cgo.SAUSAGE)
        obj.extend(center)
        obj.extend(start)
        obj.extend([
            radius,
            0.8, 0.8, 0.8,
        ])
        obj.extend(color2_list)

    coneend = cpv.add(end, cpv.scale(direction, 4.0*radius/cpv.length(direction)))
    obj.append(cgo.CONE)
    obj.extend(end)
    obj.extend(coneend)
    obj.extend([
        radius * 1.75,
        0.0,
    ])
    obj.extend(color_list * 2)
    obj.extend([
        1.0, 1.0, # Caps
    ])
    cmd.load_cgo(obj, get_unused_name('oriVec'), zoom=0)

def cafit_orientation(selection, state=STATE, visualize=1, guide=1, quiet=1):
    '''
DESCRIPTION

    Get the center and direction of a peptide by least squares
    linear fit on CA atoms.

USAGE

    cafit_orientation selection [, visualize ]

NOTES

    Requires python module "numpy".

SEE ALSO

    helix_orientation
    '''
    try:
        import numpy
    except ImportError:
        print(' Error: numpy not available')
        raise CmdException

    state, visualize, quiet = int(state), int(visualize), int(quiet)

    if int(guide):
        selection = '(%s) and guide' % (selection)

    coords = []
    cmd.iterate_state(state, selection,
            'coords.append([x,y,z])', space=locals())
    x = numpy.array(coords)

    center = x.mean(0).tolist()
    U,s,Vh = numpy.linalg.svd(x - center)

    vec = cpv.normalize(Vh[0])
    if cpv.dot_product(vec, x[-1] - x[0]) < 0:
        vec = cpv.negate(vec)

    _common_orientation(selection, center, vec, visualize, s[0], quiet)
    return center, vec

def loop_orientation(selection, state=STATE, visualize=1, quiet=1):
    '''
DESCRIPTION

    Get the center and approximate direction of a peptide. Works for any
    secondary structure.
    Averages direction of N(i)->C(i) pseudo bonds.

USAGE

    loop_orientation selection [, visualize ]

SEE ALSO

    helix_orientation
    '''
    state, visualize, quiet = int(state), int(visualize), int(quiet)

    coords = dict()
    cmd.iterate_state(state, '(%s) and name N+C' % (selection),
            'coords.setdefault(chain + resi, {})[name] = x,y,z', space=locals())

    vec = cpv.get_null()
    center = cpv.get_null()

    count = 0
    for x in coords.values():
        if 'C' in x and 'N' in x:
            vec = cpv.add(vec, cpv.sub(x['C'], x['N']))
        for coord in x.values():
            center = cpv.add(center, coord)
            count += 1

    if count == 0:
        print('warning: count == 0')
        raise CmdException

    vec = cpv.normalize(vec)
    center = cpv.scale(center, 1./count)

    _common_orientation(selection, center, vec, visualize, 2.0*len(coords), quiet)
    return center, vec

def helix_orientation(selection, state=STATE, visualize=1, cutoff=3.5, quiet=1):
    '''
DESCRIPTION

    Get the center and direction of a helix as vectors. Will only work
    for alpha helices and gives slightly different results than
    cafit_orientation. Averages direction of C(i)->O(i)->N(i+4).

USAGE

    helix_orientation selection [, visualize [, cutoff ]]

ARGUMENTS

    selection = string: atom selection of helix

    visualize = 0 or 1: show fitted vector as arrow {default: 1}

    cutoff = float: maximal hydrogen bond distance {default: 3.5}

SEE ALSO

    angle_between_helices, loop_orientation, cafit_orientation
    '''
    state, visualize, quiet = int(state), int(visualize), int(quiet)
    cutoff = float(cutoff)

    atoms = {'C': dict(), 'O': dict(), 'N': dict()}
    cmd.iterate_state(state, '(%s) and name N+O+C' % (selection),
            'atoms[name][resv] = x,y,z', space={'atoms': atoms})

    vec_list = []
    for resi in atoms['C']:
        resi_other = resi + 4
        try:
            aC = atoms['C'][resi]
            aO = atoms['O'][resi]
            aN = atoms['N'][resi_other]
        except KeyError:
            continue

        dist = cpv.distance(aN, aO)
        dist_weight = 1. - (2.8 - dist)
        angle = cpv.get_angle_formed_by(aC, aO, aN)
        angle_weight = 1. - (3.1 - angle)

        if dist_weight > 0.0 and angle_weight > 0.0:
            if not quiet:
                print(' weight:', angle_weight * dist_weight)
            vec = cpv.scale(cpv.sub(aN, aC), angle_weight * dist_weight)
            vec_list.append(vec)

    if len(vec_list) == 0:
        print('warning: count == 0')
        raise CmdException

    center = cpv.scale(_vec_sum(atoms['O'].values()), 1./len(atoms['O']))
    vec = _vec_sum(vec_list)
    vec = cpv.normalize(vec)

    _common_orientation(selection, center, vec, visualize, 1.5*len(vec_list), quiet)
    return center, vec

def plane_orientation(selection, state=STATE, visualize=1, guide=0, quiet=1):
    '''
DESCRIPTION

    Fit plane (for example beta-sheet). Can also be used with
    angle_between_helices (even though this does not fit helices).

    Returns center and normal vector of plane.
    '''
    try:
        import numpy
    except ImportError:
        print(' Error: numpy not available')
        raise CmdException

    state, visualize, quiet = int(state), int(visualize), int(quiet)

    if int(guide):
        selection = '(%s) and guide' % (selection)

    coords = list()
    cmd.iterate_state(state, selection,
            'coords.append([x,y,z])', space=locals())

    if len(coords) < 3:
        print('not enough guide atoms in selection')
        raise CmdException

    x = numpy.array(coords)
    U,s,Vh = numpy.linalg.svd(x - x.mean(0))

    # normal vector of plane is 3rd principle component
    vec = cpv.normalize(Vh[2])
    if cpv.dot_product(vec, x[-1] - x[0]) < 0:
        vec = cpv.negate(vec)

    center = x.mean(0).tolist()
    _common_orientation(selection, center, vec, visualize, 4.0, quiet)

    # plane visualize
    if visualize:
        from pymol import cgo

        dir1 = cpv.normalize(Vh[0])
        dir2 = cpv.normalize(Vh[1])
        sx = [max(i/4.0, 2.0) for i in s]

        obj = [ cgo.BEGIN, cgo.TRIANGLES, cgo.COLOR, 0.5, 0.5, 0.5 ]
        for vertex in [
                cpv.scale(dir1, sx[0]),
                cpv.scale(dir2, sx[1]),
                cpv.scale(dir2, -sx[1]),
                cpv.scale(dir1, -sx[0]),
                cpv.scale(dir2, -sx[1]),
                cpv.scale(dir2, sx[1]),
                ]:
            obj.append(cgo.VERTEX)
            obj.extend(cpv.add(center, vertex))
        obj.append(cgo.END)
        cmd.load_cgo(obj, get_unused_name('planeFit'))

    return center, vec

def angle_between_helices(selection1, selection2, method='helix',
        state1=STATE, state2=STATE, visualize=1, quiet=1):
    '''
DESCRIPTION

    Calculates the angle between two helices

USAGE

    angle_between_helices selection1, selection2 [, method [, visualize]]

ARGUMENTS

    selection1 = string: atom selection of first helix

    selection2 = string: atom selection of second helix

    method = string: function to calculate orientation {default: helix_orientation}

    visualize = 0 or 1: show fitted vector as arrow {default: 1}

EXAMPLE

    fetch 2x19, async=0
    select hel1, /2x19//B/23-36/
    select hel2, /2x19//B/40-54/
    angle_between_helices hel1, hel2
    angle_between_helices hel1, hel2, cafit

SEE ALSO

    helix_orientation, loop_orientation, cafit_orientation, angle_between_domains
    '''
    import math

    state1, state2 = int(state1), int(state2)
    visualize, quiet = int(visualize), int(quiet)

    try:
        orientation = globals()[methods_sc[str(method)]]
    except KeyError:
        print('no such method:', method)
        raise CmdException

    if not int(quiet):
        print(' Using method:', orientation.__name__)

    cen1, dir1 = orientation(selection1, state1, visualize, quiet=1)
    cen2, dir2 = orientation(selection2, state2, visualize, quiet=1)

    angle = cpv.get_angle(dir1, dir2)
    angle = math.degrees(angle)

    if not quiet:
        print(' Angle: %.2f deg' % (angle))

    if visualize:
        # measurement object for angle
        center = cpv.scale(cpv.add(cen1, cen2), 0.5)
        tmp = get_unused_name('_')
        for pos in [center,
                cpv.add(center, cpv.scale(dir1, 5.0)),
                cpv.add(center, cpv.scale(dir2, 5.0))]:
            cmd.pseudoatom(tmp, pos=list(pos), state=1)
        name = get_unused_name('angle')
        cmd.angle(name, *[(tmp, i) for i in [2,1,3]])
        cmd.delete(tmp)

        cmd.zoom('(%s) or (%s)' % (selection1, selection2), 2,
                state1 if state1 == state2 else 0)

    return angle

def angle_between_domains(selection1, selection2, method='align',
        state1=STATE, state2=STATE, visualize=1, quiet=1):
    '''
DESCRIPTION

    Angle by which a molecular selection would be rotated when superposing
    on a selection2.

    Do not use for measuring angle between helices, since the alignment of
    the helices might involve a rotation around the helix axis, which will
    result in a larger angle compared to the angle between helix axes.

USAGE

    angle_between_domains selection1, selection2 [, method ]

ARGUMENTS

    selection1 = string: atom selection of first helix

    selection2 = string: atom selection of second helix

    method = string: alignment command like "align" or "super" {default: align}

EXAMPLE

    fetch 3iplA 3iplB, async=0
    select domain1, resi 1-391
    select domain2, resi 392-475
    align 3iplA and domain1, 3iplB and domain1
    angle_between_domains 3iplA and domain2, 3iplB and domain2

SEE ALSO

    align, super, angle_between_helices
    '''
    import math

    try:
        import numpy
    except ImportError:
        print(' Error: numpy not available')
        raise CmdException

    state1, state2 = int(state1), int(state2)
    visualize, quiet = int(visualize), int(quiet)

    if cmd.is_string(method):
        try:
            method = cmd.keyword[method][0]
        except KeyError:
            print('no such method:', method)
            raise CmdException

    mobile_tmp = get_unused_name('_')
    cmd.create(mobile_tmp, selection1, state1, 1,  zoom=0)
    try:
        method(mobile=mobile_tmp, target=selection2, mobile_state=1,
                target_state=state2, quiet=quiet)
        mat = cmd.get_object_matrix(mobile_tmp)
    except:
        print(' Error: superposition with method "%s" failed' % (method.__name__))
        raise CmdException
    finally:
        cmd.delete(mobile_tmp)

    try:
        # Based on transformations.rotation_from_matrix
        # Copyright (c) 2006-2012, Christoph Gohlke

        R33 = [mat[i:i+3] for i in [0,4,8]]
        R33 = numpy.array(R33, float)

        # direction: unit eigenvector of R33 corresponding to eigenvalue of 1
        w, W = numpy.linalg.eig(R33.T)
        i = w.real.argmax()
        direction = W[:, i].real

        # rotation angle depending on direction
        m = direction.argmax()
        i,j,k,l = [
            [2,1,1,2],
            [0,2,0,2],
            [1,0,0,1]][m]
        cosa = (R33.trace() - 1.0) / 2.0
        sina = (R33[i, j] + (cosa - 1.0) * direction[k] * direction[l]) / direction[m]

        angle = math.atan2(sina, cosa)
        angle = abs(math.degrees(angle))
    except:
        print(' Error: rotation from matrix failed')
        raise CmdException

    if not quiet:
        try:
            # make this import optional to support running this script standalone
            from .querying import centerofmass, gyradius
        except (ValueError, ImportError):
            gyradius = None
            try:
                # PyMOL 1.7.1.6+
                centerofmass = cmd.centerofmass
            except AttributeError:
                centerofmass = lambda s: cpv.scale(cpv.add(*cmd.get_extent(s)), 0.5)

        center1 = centerofmass(selection1)
        center2 = centerofmass(selection2)
        print(' Angle: %.2f deg, Displacement: %.2f angstrom' % (angle, cpv.distance(center1, center2)))

        if visualize:
            center1 = numpy.array(center1, float)
            center2 = numpy.array(center2, float)
            center = (center1 + center2) / 2.0

            if gyradius is not None:
                rg = numpy.array(gyradius(selection1), float)
            else:
                rg = 10.0

            h1 = numpy.cross(center2 - center1, direction)
            h2 = numpy.dot(R33, h1)
            h1 *= rg / cpv.length(h1)
            h2 *= rg / cpv.length(h2)

            for pos in [center1, center2, center1 + h1, center1 + h2]:
                cmd.pseudoatom(mobile_tmp, pos=list(pos), state=1)

            # measurement object for angle and displacement
            name = get_unused_name('measurement')
            cmd.distance(name, *['%s`%d' % (mobile_tmp, i) for i in [1,2]])
            cmd.angle(name, *['%s`%d' % (mobile_tmp, i) for i in [3,1,4]])

            # CGO arrow for axis of rotation
            visualize_orientation(direction, center1, rg, color='blue')

            cmd.delete(mobile_tmp)

    return angle

# methods for auto-completion

methods = [
    helix_orientation,
    loop_orientation,
    cafit_orientation,
    plane_orientation,
]

methods_sc = cmd.Shortcut([func.__name__ for func in methods])

# commands and tab-completion of arguments

for func in methods:
    cmd.extend(func.__name__, func)
    cmd.auto_arg[0][func.__name__] = cmd.auto_arg[0]['zoom']

for func in [angle_between_helices, angle_between_domains]:
    cmd.extend(func.__name__, func)
    cmd.auto_arg[0][func.__name__] = cmd.auto_arg[0]['align']
    cmd.auto_arg[1][func.__name__] = cmd.auto_arg[0]['align']

cmd.auto_arg[2].update([
    ('angle_between_helices', [ methods_sc, 'method', '' ]),
])

try:
    # make this import optional to support running this script standalone
    from .fitting import align_methods_sc
    cmd.auto_arg[2].update([
        ('angle_between_domains', [ align_methods_sc, 'alignment method', '' ]),
    ])
except (ValueError, ImportError):
    pass

# vi: expandtab:smarttab
