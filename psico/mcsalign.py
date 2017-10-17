'''
(c) 2016 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException

def mcsalign(mobile, target,
        mobile_state=-1, target_state=-1,
        cycles=5, timeout=10, method='', exact=0, quiet=1):
    '''
DESCRIPTION

    Align two (ligand) selections based on Maximum-Common-Substructure.

    Requires: (rdkit | indigo), csb

ARGUMENTS

    mobile = str: atom selection of mobile object

    target = str: atom selection of target object

    mobile_state = int: object state of mobile selection
    {default: -1 = current state}

    target_state = int: object state of target selection
    {default: -1 = current state}

    cycles = int: number of weight-refinement iterations for
    weighted RMS fitting {default: 5}

    timeout = int: MCS search timeout {default: 10}

    method = indigo or rdkit {default: check availability}

    exact = 0/1: match elements and bond orders {default: 0}

EXAMPLE

    fetch 3zcf 4n8t, async=0
    mcsalign /3zcf//A/HEC, /4n8t//A/HEM
    zoom /4n8t//A/HEM, animate=2, buffer=3
    '''
    from numpy import identity, dot, take
    from csb.bio.utils import distance_sq, wfit, fit

    # moving object
    m_objects = cmd.get_object_list(mobile)
    if len(m_objects) != 1:
        # If selection covers multiple objects, call "mcsalign" for every object
        for m_object in m_objects:
            mcsalign('(%s) & model %s' % (mobile, m_object), target,
                    mobile_state, target_state, cycles, timeout, method, quiet)
        return

    # get molecules from selections
    m_sdf = get_molstr(mobile, mobile_state)
    t_sdf = get_molstr(target, target_state)

    # find maximum common substructure
    m_indices, t_indices = get_mcs_indices(method, quiet, m_sdf, t_sdf, timeout, int(exact))

    if len(m_indices) < 3:
        raise CmdException('not enough atoms in MCS')

    if not int(quiet):
        print(' MCS-Align: found MCS with %d atoms (%s)' % (len(m_indices), m_objects[0]))

    # coordinates
    Y = take(cmd.get_coords(mobile, mobile_state), m_indices, 0)
    X = take(cmd.get_coords(target, target_state), t_indices, 0)

    # weighted RMS fitting
    R, t = fit(X, Y)
    for _ in range(int(cycles)):
        data = distance_sq(Y, dot(X - t, R))
        scales = 1.0 / data.clip(1e-3)
        R, t = wfit(X, Y, scales)

    # superpose
    m = identity(4)
    m[0:3,0:3] = R
    m[0:3,3] = t
    cmd.transform_object(m_objects[0], list(m.flat), mobile_state)

def get_molstr(sele, state):
    '''Export the given selection to a molfile string'''
    from chempy import io
    model = cmd.get_model(sele, state)
    for a in model.atom:
        a.symbol = a.symbol.capitalize()
    return ''.join(io.mol.toList(model))

def get_mcs_indices(method, quiet, *args, **kwargs):
    '''Find the MCS between two molecules and return the substructure
    atom indices for both molecules.
    
    @param method: empty string, rdkit, or indigo
    @param quiet: 0 or 1 (verbosity flag)
    @param m_sdf: mobile molecule as MOL string
    @param t_sdf: target molecule as MOL string
    @param timeout: timeout in seconds (only used with rdkit)

    @return: two sequences of integers with atom indices
    @rtype: tuple
    '''

    methods = {
        'rdkit': get_mcs_indices_rdkit,
        'indigo': get_mcs_indices_indigo,
    }

    if not method:
        for method in methods:
            try:
                __import__(method)
                break
            except ImportError:
                continue
        else:
            raise CmdException('neither "rdkit" nor "indigo" available')

        if not int(quiet):
            print(" Note: using method '%s'" % (method,))

    return methods[method](*args, **kwargs)

def get_mcs_indices_rdkit(m_sdf, t_sdf, timeout, exact=False):
    try:
        from rdkit.Chem.rdFMCS import FindMCS
    except ImportError:
        # backwards compatibility
        from rdkit.Chem.MCS import FindMCS
        kwargs = {'bondCompare': 'any', 'atomCompare': 'any'}
    else:
        from rdkit.Chem import rdFMCS
        kwargs = {'bondCompare': rdFMCS.BondCompare.CompareAny,
                  'atomCompare': rdFMCS.AtomCompare.CompareAny}

    if exact:
        kwargs.clear()

    from rdkit import Chem

    m_mol = Chem.MolFromMolBlock(m_sdf, False, False)
    t_mol = Chem.MolFromMolBlock(t_sdf, False, False)

    # find maximum common substructure
    timeout = int(timeout) if timeout else None
    mcs = FindMCS([m_mol, t_mol], timeout=timeout, **kwargs)

    # backwards compatibility
    if hasattr(mcs, 'completed'):
        mcs.canceled = not mcs.completed
        mcs.smartsString = mcs.smarts

    if mcs.canceled:
        raise CmdException('MCS search has timed out')

    # atom indices of match
    patt  = Chem.MolFromSmarts(mcs.smartsString)
    m_indices = m_mol.GetSubstructMatch(patt)
    t_indices = t_mol.GetSubstructMatch(patt)

    return m_indices, t_indices

def get_mcs_indices_indigo(m_sdf, t_sdf, timeout=None, exact=False):
    import indigo
    indigo = indigo.Indigo()

    m_mol = indigo.loadMolecule(m_sdf)
    t_mol = indigo.loadMolecule(t_sdf)

    # find common substructure
    arr = indigo.createArray()
    arr.arrayAdd(m_mol)
    arr.arrayAdd(t_mol)
    mcs = indigo.extractCommonScaffold(arr, 'exact' if exact else 'approx')

    # match to scaffold
    query = indigo.loadQueryMolecule(mcs.smiles())
    m_match = indigo.substructureMatcher(m_mol).match(query)
    t_match = indigo.substructureMatcher(t_mol).match(query)

    # atom indices of match
    m_atoms = [m_match.mapAtom(a) for a in query.iterateAtoms()]
    t_atoms = [t_match.mapAtom(a) for a in query.iterateAtoms()]
    m_indices = [a.index() for a in m_atoms if a is not None]
    t_indices = [a.index() for a in t_atoms if a is not None]

    return m_indices, t_indices

# pymol commands
cmd.extend('mcsalign', mcsalign)

# auto-completion
cmd.auto_arg[0]['mcsalign'] = cmd.auto_arg[0]['align']
cmd.auto_arg[1]['mcsalign'] = cmd.auto_arg[1]['align']

# vi: ts=4:sw=4:smarttab:expandtab
