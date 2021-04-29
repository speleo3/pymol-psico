'''
(c) 2011 Thomas Holder, MPI for Developmental Biology
(c) 2011 Tsjerk Wassenaar (gyradius code)

License: BSD-2-Clause
'''

from pymol import cmd, CmdException
from pymol import selector

def centerofmass(selection='(all)', state=-1, quiet=1):
    '''
DESCRIPTION

    Calculates the center of mass. Considers atom mass and occupancy.

ARGUMENTS

    selection = string: atom selection {default: all}

    state = integer: object state, -1 for current state, 0 for all states
    {default: -1}

EXAMPLE

    from psico.querying import *
    x = centerofmass('chain A')
    r = gyradius('chain A')
    cmd.pseudoatom('com', pos=x, vdw=r)

SEE ALSO

    gyradius
    '''
    from chempy import cpv
    state, quiet = int(state), int(quiet)
    if state < 0:
        states = [cmd.get_state()]
    elif state == 0:
        states = list(range(1, cmd.count_states(selection)+1))
    else:
        states = [state]
    com = cpv.get_null()
    totmass = 0.0
    for state in states:
        model = cmd.get_model(selection, state)
        for a in model.atom:
            if a.q == 0.0:
                continue
            m = a.get_mass() * a.q
            com = cpv.add(com, cpv.scale(a.coord, m))
            totmass += m
    com = cpv.scale(com, 1./totmass)
    if not quiet:
        print(' Center of Mass: [%8.3f,%8.3f,%8.3f]' % tuple(com))
    return com

def gyradius(selection='(all)', state=-1, quiet=1):
    '''
DESCRIPTION
 
    Radius of gyration

    Based on: http://pymolwiki.org/index.php/Radius_of_gyration

SEE ALSO

    centerofmass
    '''
    from chempy import cpv
    state, quiet = int(state), int(quiet)
    if state < 0:
        states = [cmd.get_state()]
    elif state == 0:
        states = list(range(1, cmd.count_states(selection)+1))
    else:
        states = [state]
    rg_sq_list = []
    for state in states:
        model = cmd.get_model(selection, state)
        x = [i.coord for i in model.atom]
        mass = [i.get_mass() * i.q for i in model.atom if i.q > 0]
        xm = [cpv.scale(v,m) for v,m in zip(x,mass)]
        tmass = sum(mass)
        rr = sum(cpv.dot_product(v,vm) for v,vm in zip(x,xm))
        mm = sum((sum(i)/tmass)**2 for i in zip(*xm))
        rg_sq_list.append(rr/tmass - mm)
    rg = (sum(rg_sq_list)/len(rg_sq_list))**0.5
    if not quiet:
        print(' Radius of gyration: %.2f' % (rg))
    return rg

def get_alignment_coords(name, active_only=0, state=-1, quiet=0):
    '''
DESCRIPTION

    API only function. Returns a dictionary with items

        (object name, Nx3 coords list)
 
    N is the number of alignment columns without gaps.
 
EXAMPLE

    import numpy
    from psico.multistuff import *
    from psico.querying import *

    extra_fit('name CA', cycles=0, object='aln')
    x = get_alignment_coords('aln')
    m = numpy.array(x.values())
    '''
    active_only, state, quiet = int(active_only), int(state), int(quiet)
    aln = cmd.get_raw_alignment(name, active_only)
    object_list = cmd.get_object_list(name)
    idx2coords = dict()
    cmd.iterate_state(state, name, 'idx2coords[model,index] = (x,y,z)',
            space={'idx2coords': idx2coords})
    allcoords = dict((model, []) for model in object_list)
    for pos in aln:
        if len(pos) != len(object_list):
            continue
        for model,index in pos:
            allcoords[model].append(idx2coords[model,index])
    return allcoords

def get_sasa(selection, state=-1, dot_density=5, quiet=1):
    '''
DESCRIPTION

    Get solvent accesible surface area

SEE ALSO

    get_area
    pymol.util.get_sasa (considered broken!)
    '''
    state, dot_density, quiet = int(state), int(dot_density), int(quiet)
    if state < 1:
        state = cmd.get_state()
    n = cmd.get_unused_name('_')
    cmd.create(n, selection, state, 1, zoom=0, quiet=1)
    cmd.set('dot_solvent', 1, n)
    if dot_density > -1:
        cmd.set('dot_density', dot_density, n)
    r = cmd.get_area(n, quiet=int(quiet))
    cmd.delete(n)
    return r

def get_sasa_ball(selection, state=-1, quiet=1):
    '''
DESCRIPTION

    Get solvent accesible surface area using BALL.NumericalSAS

    http://www.ball-project.org/
    '''
    import BALL
    import tempfile, os

    state, quiet = int(state), int(quiet)
    radius = cmd.get_setting_float('solvent_radius')

    filename = tempfile.mktemp('.pdb')
    cmd.save(filename, selection, state, 'pdb')
    system = BALL.System()
    BALL.PDBFile(filename) >> system
    os.remove(filename)

    fragment_db = BALL.FragmentDB('')
    system.apply(fragment_db.normalize_names)
    system.apply(BALL.AssignRadiusProcessor('radii/PARSE.siz'))

    sas = BALL.NumericalSAS()
    sas_options = BALL.Options()
    sas_options.setBool(sas.Option.COMPUTE_AREA, True)
    sas_options.setBool(sas.Option.COMPUTE_SURFACE, False)
    sas_options.setReal(sas.Option.PROBE_RADIUS, radius)
    sas.setOptions(sas_options)
    sas(system)
    area = sas.getTotalArea()
    if not quiet:
        print(' get_sasa_ball: %.3f Angstroms^2.' % (area))
    return area

def get_sasa_mmtk(selection, state=-1, hydrogens='auto', quiet=1):
    '''
DESCRIPTION

    Get solvent accesible surface area using MMTK.MolecularSurface

    http://dirac.cnrs-orleans.fr/MMTK/

    This command is very picky with missing atoms and wrong atom naming.

SEE ALSO

    stub2ala, get_sasa, get_sasa_ball
    '''
    import MMTK
    from MMTK.PDB import PDBConfiguration
    from MMTK.Proteins import Protein
    from MMTK.MolecularSurface import surfaceAndVolume

    try:
        from cStringIO import StringIO
    except ImportError:
        from io import StringIO

    selection = selector.process(selection)
    state, quiet = int(state), int(quiet)
    radius = cmd.get_setting_float('solvent_radius')

    if hydrogens == 'auto':
        if cmd.count_atoms('(%s) and hydro' % selection) > 0:
            hydrogens = 'all'
        else:
            hydrogens = 'no_hydrogens'
    elif hydrogens == 'none':
        hydrogens = 'no_hydrogens'

    conf = PDBConfiguration(StringIO(cmd.get_pdbstr(selection)))
    system = Protein(conf.createPeptideChains(hydrogens))

    try:
        area, volume = surfaceAndVolume(system, radius * 0.1)
    except:
        print(' Error: MMTK.MolecularSurface.surfaceAndVolume failed')
        raise CmdException

    if not quiet:
        print(' get_sasa_mmtk: %.3f Angstroms^2 (volume: %.3f Angstroms^3).' % (area * 1e2, volume * 1e3))
    return area * 1e2

def get_raw_distances(names='', state=1, selection='all', quiet=1):
    '''
DESCRIPTION

    Get the list of pair items from distance objects. Each list item is a
    tuple of (index1, index2, distance).

    Based on a script from Takanori Nakane, posted on pymol-users mailing list.
    http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10143.html

ARGUMENTS

    names = string: names of distance objects (no wildcards!) {default: all
    measurement objects}

    state = integer: object state {default: 1}

    selection = string: atom selection {default: all}

SEE ALSO

    select_distances, cmd.find_pairs, cmd.get_raw_alignment
    '''
    from chempy import cpv

    state, quiet = int(state), int(quiet)
    if state < 1:
        state = cmd.get_state()

    valid_names = cmd.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                raise CmdException('no such distance object: ' + name)

    raw_objects = cmd.get_session(names, 1, 1, 0, 0)['names']

    xyz2idx = {}
    cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,index)',
            space=locals())

    r = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state-1][1]
            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i+3])
            xyz2 = tuple(points[i+3:i+6])
            try:
                r.append((xyz2idx[xyz1], xyz2idx[xyz2], cpv.distance(xyz1, xyz2)))
                if not quiet:
                    print(' get_raw_distances: ' + str(r[-1]))
            except KeyError:
                if quiet < 0:
                    print(' Debug: no index for %s %s' % (xyz1, xyz2))
    return r

def get_color(selection, which=0, mode=0):
    '''
DESCRIPTION

    API only. Returns the color of the first/middle/... guide atom in
    selection.

ARGUMENTS

    which = 0: color of first atom
    which = 1: color of middle atom
    which = 2: most frequent color

    mode = 0: color index or color string
    mode = 1: color tuple
    mode = 2: color string in hash-hex format (for HTML, matplotlib, ...)
    '''
    s_first = 'first' if which == 0 else ''

    try:
        colors = []

        for s_guide in ('guide', 'elem C', 'all'):
            cmd.iterate('{} (({}) & {})'.format(s_first, selection, s_guide),
                    'colors.append(color)', space=locals())
            if colors:
                break

        if which == 2:
            color = max((colors.count(color), color) for color in colors)[1]
        else:
            color = colors[len(colors) // 2]

        if color >= 0x40000000:
            color = '0x%06x' % (color & 0xFFFFFF)
    except:
        print(' Warning: could not get color for ' + str(selection))
        color = 'gray'
    if mode > 0:
        color = cmd.get_color_tuple(color)
    if mode == 2:
        return '#%02x%02x%02x' % tuple(int(0xFF * v) for v in color)
    return color

def get_object_name(selection, strict=0):
    '''
DESCRIPTION

    Returns the object name for given selection.
    '''
    names = cmd.get_object_list('(' + selection + ')')
    if len(names) == 0:
        raise CmdException('No objects in selection')
    if strict and len(names) > 1:
        raise CmdException('Selection spans more than one object')
    return names[0]

def get_object_state(name):
    '''
DESCRIPTION

    Returns the effective object state.
    '''
    states = cmd.count_states(name)
    if states < 2 and cmd.get_setting_boolean('static_singletons'):
        return 1
    state = cmd.get_setting_int('state', name)
    if state > states:
        raise CmdException('Invalid state %d for object %s' % (state, name))
    return state

def get_selection_state(selection):
    '''
DESCRIPTION

    Returns the effective object state for all objects in given selection.
    Raises exception if objects are in different states.
    '''
    state_set = set(map(get_object_state,
        cmd.get_object_list('(' + selection + ')')))
    if len(state_set) != 1:
        if len(state_set) == 0:
            return 1
        raise CmdException('Selection spans multiple object states')
    return state_set.pop()

def _get_coords(selection, state=-1):
    '''
DESCRIPTION

    API only. Returns the (natoms, 3) coordinate matrix for a given state.
    Considers the object rotation matrix. 
    '''
    if state < 0:
        state = get_selection_state(selection)
    return cmd.get_model(selection, state).get_coord_list()

try:
    # new in PyMOL 1.7.4
    from pymol.querying import get_coords
except ImportError:
    # fallback
    get_coords = _get_coords

def get_ensemble_coords(selection):
    '''
DESCRIPTION

    API only. Returns the (nstates, natoms, 3) coordinate matrix. Considers
    the object rotation matrix. 
    '''
    nstates = cmd.count_states(selection)
    if get_coords is not _get_coords:
        # PyMOL 1.7.4 implementation
        return get_coords(selection, 0).reshape((nstates, -1, 3))
    return [get_coords(selection, state)
            for state in range(1, nstates + 1)]

if 'centerofmass' not in cmd.keyword:
    cmd.extend('centerofmass', centerofmass)
cmd.extend('gyradius', gyradius)
cmd.extend('get_sasa', get_sasa)
cmd.extend('get_sasa_ball', get_sasa_ball)
cmd.extend('get_sasa_mmtk', get_sasa_mmtk)
cmd.extend('get_raw_distances', get_raw_distances)

cmd.auto_arg[0].update([
    ('centerofmass', cmd.auto_arg[0]['zoom']),
    ('gyradius', cmd.auto_arg[0]['zoom']),
    ('get_sasa', cmd.auto_arg[0]['zoom']),
    ('get_sasa_ball', cmd.auto_arg[0]['zoom']),
    ('get_sasa_mmtk', cmd.auto_arg[0]['zoom']),
    ('get_raw_distances', [
        lambda: cmd.Shortcut(cmd.get_names_of_type('object:measurement')),
        'distance object', '']),
])

# vi: ts=4:sw=4:smarttab:expandtab
