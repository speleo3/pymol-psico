'''
(c) 2011 Thomas Holder, MPI for Developmental Biology
(c) 2011 Tsjerk Wassenaar (gyradius code)

License: BSD-2-Clause
'''

from pymol import cmd, CmdException
from pymol import selector

CURRENT_STATE = -1  # pymol.constants.CURRENT_STATE


def centerofmass(selection='(all)', state=-1, quiet=1, *, _self=cmd):
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
        states = [_self.get_state()]
    elif state == 0:
        states = list(range(1, _self.count_states(selection) + 1))
    else:
        states = [state]
    com = cpv.get_null()
    totmass = 0.0
    for state in states:
        model = _self.get_model(selection, state)
        for a in model.atom:
            if a.q == 0.0:
                continue
            m = a.get_mass() * a.q
            com = cpv.add(com, cpv.scale(a.coord, m))
            totmass += m
    com = cpv.scale(com, 1. / totmass)
    if not quiet:
        print(' Center of Mass: [%8.3f,%8.3f,%8.3f]' % tuple(com))
    return com


def gyradius(selection='(all)', state=-1, quiet=1, *, _self=cmd):
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
        states = [_self.get_state()]
    elif state == 0:
        states = list(range(1, _self.count_states(selection) + 1))
    else:
        states = [state]
    rg_sq_list = []
    for state in states:
        model = _self.get_model(selection, state)
        x = [i.coord for i in model.atom]
        mass = [i.get_mass() * i.q for i in model.atom if i.q > 0]
        xm = [cpv.scale(v, m) for v, m in zip(x, mass)]
        tmass = sum(mass)
        rr = sum(cpv.dot_product(v, vm) for v, vm in zip(x, xm))
        mm = sum((sum(i) / tmass)**2 for i in zip(*xm))
        rg_sq_list.append(rr / tmass - mm)
    rg = (sum(rg_sq_list) / len(rg_sq_list))**0.5
    if not quiet:
        print(' Radius of gyration: %.2f' % (rg))
    return rg


def get_alignment_coords(name, active_only=0, state=-1, quiet=0, *, _self=cmd):
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
    aln = _self.get_raw_alignment(name, active_only)
    object_list = _self.get_object_list(name)
    idx2coords = dict()
    _self.iterate_state(state, name, 'idx2coords[model,index] = (x,y,z)',
            space={'idx2coords': idx2coords})
    allcoords = dict((model, []) for model in object_list)
    for pos in aln:
        if len(pos) != len(object_list):
            continue
        for model, index in pos:
            allcoords[model].append(idx2coords[model, index])
    return allcoords


def get_sasa(selection, state=-1, dot_density=5, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Get solvent accesible surface area

SEE ALSO

    get_area
    pymol.util.get_sasa (considered broken!)
    '''
    state, dot_density, quiet = int(state), int(dot_density), int(quiet)
    if state < 1:
        state = _self.get_state()
    n = _self.get_unused_name('_')
    _self.create(n, selection, state, 1, zoom=0, quiet=1)
    _self.set('dot_solvent', 1, n)
    if dot_density > -1:
        _self.set('dot_density', dot_density, n)
    r = _self.get_area(n, quiet=int(quiet))
    _self.delete(n)
    return r


def get_sasa_ball(selection, state=-1, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Get solvent accesible surface area using BALL.NumericalSAS

    http://www.ball-project.org/
    '''
    import BALL
    import tempfile, os

    state, quiet = int(state), int(quiet)
    radius = _self.get_setting_float('solvent_radius')

    filename = tempfile.mktemp('.pdb')
    _self.save(filename, selection, state, 'pdb')
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


def get_sasa_mmtk(selection, state=-1, hydrogens='auto', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Get solvent accesible surface area using MMTK.MolecularSurface

    http://dirac.cnrs-orleans.fr/MMTK/

    This command is very picky with missing atoms and wrong atom naming.

SEE ALSO

    stub2ala, get_sasa, get_sasa_ball
    '''
    from MMTK.PDB import PDBConfiguration
    from MMTK.Proteins import Protein
    from MMTK.MolecularSurface import surfaceAndVolume
    from io import StringIO

    selection = selector.process(selection)
    state, quiet = int(state), int(quiet)
    radius = _self.get_setting_float('solvent_radius')

    if hydrogens == 'auto':
        if _self.count_atoms('(%s) and hydro' % selection) > 0:
            hydrogens = 'all'
        else:
            hydrogens = 'no_hydrogens'
    elif hydrogens == 'none':
        hydrogens = 'no_hydrogens'

    conf = PDBConfiguration(StringIO(_self.get_pdbstr(selection)))
    system = Protein(conf.createPeptideChains(hydrogens))

    try:
        area, volume = surfaceAndVolume(system, radius * 0.1)
    except:
        print(' Error: MMTK.MolecularSurface.surfaceAndVolume failed')
        raise CmdException

    if not quiet:
        print(' get_sasa_mmtk: %.3f Angstroms^2 (volume: %.3f Angstroms^3).' % (area * 1e2, volume * 1e3))
    return area * 1e2


def get_raw_distances(names='', state=1, selection='all', quiet=1, *, _self=cmd):
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
        state = _self.get_state()

    valid_names = _self.get_names_of_type('object:measurement')
    if names == '':
        names = ' '.join(valid_names)
    else:
        for name in names.split():
            if name not in valid_names:
                raise CmdException('no such distance object: ' + name)

    raw_objects = _self.get_session(names, 1, 1, 0, 0)['names']

    xyz2idx = {}
    _self.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,index)',
            space=locals())

    r = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state - 1][1]
            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                r.append((xyz2idx[xyz1], xyz2idx[xyz2], cpv.distance(xyz1, xyz2)))
                if not quiet:
                    print(' get_raw_distances: ' + str(r[-1]))
            except KeyError:
                if quiet < 0:
                    print(' Debug: no index for %s %s' % (xyz1, xyz2))
    return r


def get_color(selection, which=0, mode=0, *, _self=cmd):
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
            _self.iterate('{} (({}) & {})'.format(s_first, selection, s_guide),
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
        color = _self.get_color_tuple(color)
    if mode == 2:
        return '#%02x%02x%02x' % tuple(int(0xFF * v) for v in color)
    return color


def get_object_name(selection, strict=0, *, _self=cmd):
    '''
DESCRIPTION

    Returns the object name for given selection.
    '''
    names = _self.get_object_list('(' + selection + ')')
    if len(names) == 0:
        raise CmdException('No objects in selection')
    if strict and len(names) > 1:
        raise CmdException('Selection spans more than one object')
    return names[0]


def get_object_state(name, *, _self=cmd):
    '''
DESCRIPTION

    Returns the effective object state.
    '''
    states = _self.count_states(name)
    if states < 2 and _self.get_setting_boolean('static_singletons'):
        return 1
    state = _self.get_setting_int('state', name)
    if state > states:
        raise CmdException('Invalid state %d for object %s' % (state, name))
    return state


def get_selection_state(selection, *, _self=cmd):
    '''
DESCRIPTION

    Returns the effective object state for all objects in given selection.
    Raises exception if objects are in different states.
    '''
    state_set = set(map(get_object_state, _self.get_object_list(selection)))
    if len(state_set) != 1:
        if len(state_set) == 0:
            return 1
        raise CmdException('Selection spans multiple object states')
    return state_set.pop()


def get_ensemble_coords(selection, *, _self=cmd):
    '''
DESCRIPTION

    API only. Returns the (nstates, natoms, 3) coordinate matrix. Considers
    the object rotation matrix.
    '''
    nstates = _self.count_states(selection)
    return _self.get_coords(selection, 0).reshape((nstates, -1, 3))


def iterate_to_list(selection, expression, *, space=None, _self=cmd):
    """
    API-only function to capture "iterate" results in a list.
    """
    outlist = []
    _self.iterate(selection, "outlist.append(({}))".format(expression),
            space=dict(space or (), outlist=outlist))
    return outlist


def iterate_state_to_list(state: int,
                          selection: str,
                          expression: str,
                          *,
                          space=None,
                          _self=cmd) -> list:
    """
    API-only function to capture "iterate_state" results in a list.
    """
    outlist = []
    _self.iterate_state(state,
                        selection,
                        f"outlist.append(({expression}))",
                        space=dict(space or (), outlist=outlist))
    return outlist


def iterate_to_sele(selection: str,
                    expression: str,
                    *,
                    space=None,
                    _self=cmd) -> str:
    """
    API-only function to get a selection expression for "iterate" results which
    evaluate to True.
    """
    space = dict(space or ())
    space["ids"] = []

    _self.iterate(selection,
                  f"ids.append((model,index)) if ({expression}) else None",
                  space=space)

    return " ".join(f"{model}`{index}" for (model, index) in space["ids"])


def csp(sele1, sele2='', quiet=1, var="formal_charge", _self=cmd):
    """
DESCRIPTION

    Charge Symmetry Parameter between two selections. Can be used to compute
    FvCSP according to Sharma 2014.

    If only sele1 is given, it must contain excatly two chains.
    """
    if not sele2:
        chains = _self.get_chains(sele1)

        if len(chains) != 2:
            raise CmdException("need two chains")

        sele2 = '({}) & chain "{}"'.format(sele1, chains[1])
        sele1 = '({}) & chain "{}"'.format(sele1, chains[0])

    charges1 = iterate_to_list(sele1, var, _self=_self)
    charges2 = iterate_to_list(sele2, var, _self=_self)

    r = sum(charges1) * sum(charges2)

    if not int(quiet):
        print(" csp: {}".format(r))

    return r


def extinction_coefficient(selection="all", state=-1, *, quiet=1, _self=cmd):
    """
DESCRIPTION

    Extinction coefficient at 280 nm.
    """
    from pymol.util import compute_mass

    nW = _self.count_atoms(f"({selection}) & resn TRP & guide", state=state)
    nY = _self.count_atoms(f"({selection}) & resn TYR & guide", state=state)
    nSS = _self.count_atoms(
        f"({selection}) & resn CYS & elem S & bound_to elem S",
        state=state) // 2
    eps = nW * 5500 + nY * 1490 + nSS * 125
    implicit = _self.count_atoms(f"({selection}) & hydro") == 0
    mass = compute_mass(selection, state, implicit=implicit, _self=_self)
    A_280 = eps / mass

    if not int(quiet):
        print(" Extinction coefficient at 280nm: "
              f"{eps}/(M*cm), {A_280:.4f} g/L")

    return (eps, A_280)


@cmd.extendaa(cmd.auto_arg[1]['distance'], cmd.auto_arg[1]['distance'])
def shortest_distance(selection1: str,
                      selection2: str,
                      state1: int = CURRENT_STATE,
                      state2: int = CURRENT_STATE,
                      name: str = "shortest",
                      *,
                      quiet: int = 1,
                      _self=cmd):
    '''
DESCRIPTION

    Finds the shortest pairwise distance between two selections.

ARGUMENTS

    selection1 = string: first atom selection

    selection2 = string: second atom selection

    state1 = state of selection1 {default: current state}

    state2 = state of selection2 {default: current state}

    name = string: name of the object to create {default: shortest}

    quiet = 0 or 1: print results to the terminal {default: 1}

EXAMPLE

    fetch 2xwu
    shortest_distance chain A, chain B
    '''
    from math import sqrt
    from chempy import cpv

    class ShortestDistanceAtom:

        def __init__(self, model: str, segi: str, chain: str, resn: str,
                     resi: int, name: str, index: int, coord, state: int) -> None:
            self.model = model
            self.segi = segi
            self.chain = chain
            self.resn = resn
            self.resi = resi
            self.name = name
            self.index = index
            self.coord = coord
            self.state = state

        def asSelection(self) -> str:
            return f"/{self.model}/{self.segi}/{self.chain}/{self.resn}`{self.resi}/{self.name}"

        def __eq__(self, other) -> bool:
            if isinstance(other, ShortestDistanceAtom):
                return (self.model, self.index, self.state) == \
                       (other.model, other.index, other.state)
            return False

    def get_atoms(state: int, sele: str):
        return iterate_state_to_list(state, sele,
            "ShortestDistanceAtom(model, segi, chain, resn, resi, name, index, (x, y, z), state)",
            space={'ShortestDistanceAtom': ShortestDistanceAtom}, _self=_self)

    sele_1_atoms = get_atoms(state1, selection1)
    sele_2_atoms = get_atoms(state2, selection2)

    # Calculate the shortest distance
    min_distance_sq = None
    closest_pair = None
    for sele_a_atom in sele_1_atoms:
        for sele_b_atom in sele_2_atoms:
            if sele_a_atom == sele_b_atom:
                continue
            dist_sq = cpv.distance_sq(sele_a_atom.coord, sele_b_atom.coord)
            if min_distance_sq is None or dist_sq < min_distance_sq:
                min_distance_sq = dist_sq
                closest_pair = (sele_a_atom, sele_b_atom)

    if closest_pair is None:
        raise ValueError("No atoms found in selections")

    # Show the shortest distance as a dashed line
    sele_1 = closest_pair[0].asSelection()
    sele_2 = closest_pair[1].asSelection()
    name = _self.get_unused_name(name)
    _self.distance(name, sele_1, sele_2, quiet=quiet)

    min_distance = sqrt(min_distance_sq)
    if not quiet:
        print(
            f"Shortest distance: {min_distance:.2f} Å between {sele_1} and {sele_2}")
    return (min_distance, sele_1, sele_2)


@cmd.extend
def isoelectric_point(selection: str = "polymer",
                      *,
                      ph: float = 7,
                      quiet: bool = 1,
                      _self=cmd):
    """
DESCRIPTION

    Compute isoelectric point and charge at given pH.
    """
    import io
    from Bio import SeqIO
    from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

    ph, quiet = float(ph), int(quiet)
    label = ""

    fasta = _self.get_fastastr(selection)
    if not fasta:
        if not quiet:
            print(" empty selection")
        return

    records = list(SeqIO.parse(io.StringIO(fasta), "fasta"))

    if len(records) > 1 and not quiet:
        label = " combined:"
        for rec in records:
            protein = IsoelectricPoint(rec.seq)
            print(f" {rec.id}: pI={protein.pi():.2f}"
                  f" charge={protein.charge_at_pH(ph):.2f}")

    seq = "".join(str(rec.seq) for rec in records)
    protein = IsoelectricPoint(seq)
    pi = protein.pi()

    if not quiet:
        print(f"{label} pI={pi:.2f} charge={protein.charge_at_pH(ph):.2f}")

    return pi


@cmd.extend
def get_segis(selection="all", *, quiet=1, _self=cmd) -> set:
    """
DESCRIPTION

    Get the set of segment identifiers.
    """
    segis = set(iterate_to_list(selection, "segi", _self=_self))

    if not int(quiet):
        print(f" get_segis: {segis}")

    return segis


if 'centerofmass' not in cmd.keyword:
    cmd.extend('centerofmass', centerofmass)
cmd.extend('gyradius', gyradius)
cmd.extend('get_sasa', get_sasa)
cmd.extend('get_sasa_ball', get_sasa_ball)
cmd.extend('get_sasa_mmtk', get_sasa_mmtk)
cmd.extend('get_raw_distances', get_raw_distances)
cmd.extend('csp', csp)
cmd.extend('extinction_coefficient', extinction_coefficient)

cmd.auto_arg[0].update([
    ('centerofmass', cmd.auto_arg[0]['zoom']),
    ('gyradius', cmd.auto_arg[0]['zoom']),
    ('get_sasa', cmd.auto_arg[0]['zoom']),
    ('get_segis', cmd.auto_arg[0]['zoom']),
    ('get_sasa_ball', cmd.auto_arg[0]['zoom']),
    ('get_sasa_mmtk', cmd.auto_arg[0]['zoom']),
    ('get_raw_distances', [
        lambda: cmd.Shortcut(cmd.get_names_of_type('object:measurement')),
        'distance object', '']),
    ('csp', cmd.auto_arg[0]['zoom']),
])

# vi: ts=4:sw=4:smarttab:expandtab
