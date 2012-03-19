'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def split_chains(selection='(all)', prefix=None):
    '''
DESCRIPTION

    Create a single object for each chain in selection

SEE ALSO

    split_states
    '''
    count = 0
    models = cmd.get_object_list('(' + selection + ')')
    for model in models:
        for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
            count += 1
            if not prefix:
                name = '%s_%s' % (model, chain)
            else:
                name = '%s%04d' % (prefix, count)
            cmd.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
        cmd.disable(model)

def rmsf2b(selection='name CA', linearscale=1.0, quiet=1):
    '''
DESCRIPTION

    Determine the root mean square fluctuation (RMSF) per atom for a
    multi-state object and assign b-factor

ARGUMENTS

    selection = string: atom selection {default: name CA}

    linearscale = float: if linearscale <= 0, then use real b-factor equation,
    else use b=(rmsf*linearscale) {default: 1.0}

SEE ALSO

    spheroid, rmsf_states.py from Robert Campbell
    '''
    from numpy import array, sqrt, pi
    linearscale = float(linearscale)
    n_atoms = cmd.count_atoms(selection)
    n_states = cmd.count_states(selection)
    if n_atoms == 0 or n_states < 2:
        print ' Error: not enough atoms or states'
        raise CmdException
    coords = []
    for state in range(1, n_states + 1):
        state_coords = cmd.get_model(selection, state).get_coord_list()
        if len(state_coords) != n_atoms:
            print ' Error: number of atoms in states not equal'
            raise CmdException
        coords.append(state_coords)
    coords = array(coords)
    u_sq = coords.var(0).sum(1) # var over states, sum over x,y,z
    b_array = sqrt(u_sq) * linearscale if linearscale > 0.0 \
            else 8 * pi**2 * u_sq
    cmd.alter(selection, 'b = b_iter.next()', space={'b_iter': iter(b_array)})
    if not int(quiet):
        print ' Average RMSF: %.2f' % (sqrt(u_sq).mean())
    return b_array

def set_sequence(sequence, selection='all', start=1):
    '''
DESCRIPTION

    Alters the residue names according to given sequence

ARGUMENTS

    sequence = string: amino acid sequence in one-letter code

    selection = string: atom selection {default: all}

    start = int: residue number to start from {default: 1}
    '''
    import re
    from . import three_letter
    sequence = re.sub(r'\s+', '', sequence)
    start = int(start)
    for i, aa in enumerate(sequence):
        cmd.alter('(%s) and resi %d' % (selection, i+start),
                'resn=' + repr(three_letter.get(aa.upper(), 'UNK')))

def alphatoall(selection='polymer', properties='b', operator='byca', quiet=1):
    '''
DESCRIPTION

    Expand any given property of the CA atoms to all atoms in the residue

    Enhanced version of http://pymolwiki.org/index.php/AlphaToAll

ARGUMENTS

    selection = string: atom selection {default: polymer}

    properties = string: space separated list of atom properties {default: b}
    '''
    properties = '(' + ','.join(properties.split()) + ')'
    space = {'props': dict()}
    cmd.iterate('%s (%s)' % (operator, selection), 'props[model,segi,chain,resi] = ' + properties,
            space=space)
    cmd.alter(selection,
            properties + ' = props.get((model,segi,chain,resi), ' + properties + ')',
            space=space)
    if not int(quiet):
        print ' Modified %d residues' % (len(space['props']))

def mse2met(selection='all', quiet=1):
    '''
DESCRIPTION

    Mutate selenomethionine to methionine
    '''
    quiet = int(quiet)
    x = cmd.alter('(%s) and MSE/SE' % selection, 'name="SD";elem="S"')
    cmd.flag('ignore', '(%s) and MSE/' % (selection), 'clear')
    cmd.alter('(%s) and MSE/' % selection, 'resn="MET";type="ATOM"')
    if not quiet:
        print 'Altered %d MSE residues to MET' % (x)
    cmd.sort()

def polyala(selection='all', quiet=1):
    '''
DESCRIPTION

    Mutate any residue to Alanine (except Glycines)

SEE ALSO

    stub2ala
    '''
    quiet = int(quiet)
    cmd.remove('polymer and (%s) and not name C+N+O+CA+CB+OXT' % (selection))
    cmd.alter('polymer and (%s) and not resn GLY' % (selection), 'resn = "ALA"')
    cmd.sort()

def stub2ala(selection='all', quiet=1):
    '''
DESCRIPTION

    Mutate stub residues to ALA

SEE ALSO

    polyala
    '''
    quiet = int(quiet)
    namesets = dict()
    lookslike = {
        # keys are sorted tuples of backbone atoms
        ('CA',): 'GLY',
        ('CA', 'CB'): 'ALA',
        ('CA', 'CB', 'CG1', 'CG2'): 'VAL',
        ('CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ'): 'PHE',
    }
    cmd.iterate('(%s) and polymer and (not hydro) and (not name C+N+O+OXT)' % (selection),
            'namesets.setdefault((model,segi,chain,resv,resn,resi), set()).add(name)',
            space={'namesets': namesets, 'set': set})
    for key in sorted(namesets):
        resn = key[-2]
        name_tuple = tuple(sorted(namesets[key]))
        key_str = '/%s/%s/%s/%s`%s' % (key[:3] + key[4:])
        if name_tuple == ('CA', 'CB', 'CG'):
            key_str_cg = key_str + '/CG'
            if not quiet:
                print 'Removing', key_str_cg
            cmd.remove(key_str_cg)
            name_tuple = ('CA', 'CB')
        lookslike_resn = lookslike.get(name_tuple, resn)
        if lookslike_resn != resn:
            if not quiet:
                print 'Altering %s to %s' % (key_str, lookslike_resn)
            cmd.alter(key_str, 'resn = %s' % (repr(lookslike_resn)))
    cmd.sort()

def remove_alt(selection='all', keep='A', quiet=1):
    '''
DESCRIPTION

    Remove alternative location atoms.

USAGE

    remove_alt [selection [, keep]]

ARGUMENTS

    selection = string: atom selection

    keep = string: AltLoc to keep {default: A}
    '''
    cmd.remove('(%s) and not alt +%s' % (selection, keep), quiet=int(quiet))
    cmd.alter(selection, '(alt,q)=("",1.0)')
    cmd.sort()

def _common_ss_alter(selection, ss_dict, ss_map, raw=''):
    '''
DESCRIPTION

    Shared code of 'dssp' and 'stride' functions.
    '''
    if raw != 'ss':
        cmd.alter(selection, 'ss = ss_map.get(ss_dict.get((model,chain,resi)), "")',
                space={'ss_dict': ss_dict, 'ss_map': ss_map})
    if raw != '':
        cmd.alter(selection, raw + ' = ss_dict.get((model,chain,resi), "")',
                space={'ss_dict': ss_dict})
    cmd.rebuild(selection, 'cartoon')

def dssp(selection='(all)', exe='', raw='', state=-1, quiet=1):
    '''
DESCRIPTION

    Secondary structure assignment with DSSP.
    http://swift.cmbi.ru.nl/gv/dssp/

ARGUMENTS

    selection = string: atom selection {default: all}

    exe = string: name of dssp executable {default: dsspcmbi}

    raw = string: atom property to load raw dssp class into {default: ''}

EXAMPLE

    dssp all, /sw/bin/dsspcmbi, raw=text_type
    color gray
    color red, text_type H
    color orange, text_type G
    color yellow, text_type E
    color wheat, text_type B
    color forest, text_type T
    color green, text_type S
    set cartoon_discrete_colors, 1

SEE ALSO

    dss, stride
    '''
    from subprocess import Popen, PIPE
    import tempfile, os

    state, quiet = int(state), int(quiet)

    if exe == '':
        from . import which
        exe = which('dsspcmbi', 'dssp', 'dssp-2')
    ss_map = {
        'B': 'S', # residue in isolated beta-bridge
        'E': 'S', # extended strand, participates in beta ladder
        'T': 'L', # hydrogen bonded turn
        'G': 'H', # 3-helix (3/10 helix)
        'H': 'H', # alpha helix
        'I': 'H', # 5 helix (pi helix)
        'S': 'L', # bend
        ' ': 'L', # loop or irregular
    }
    tmpfilepdb = tempfile.mktemp('.pdb')
    ss_dict = dict()
    for model in cmd.get_object_list(selection):
        cmd.save(tmpfilepdb, '%s and (%s)' % (model, selection), state)
        try:
            process = Popen([exe, tmpfilepdb], stdout=PIPE)
        except OSError:
            print 'Error: Cannot execute exe=' + exe
            raise CmdException
        for line in process.stdout:
            if line.startswith('  #  RESIDUE'):
                break
        for line in process.stdout:
            resi = line[5:11].strip()
            chain = line[11].strip()
            ss = line[16]
            ss_dict[model,chain,resi] = ss
    os.remove(tmpfilepdb)
    _common_ss_alter(selection, ss_dict, ss_map, raw)

def stride(selection='(all)', exe='stride', raw='', state=-1, quiet=1):
    '''
DESCRIPTION

    Secondary structure assignment with STRIDE.
    http://webclu.bio.wzw.tum.de/stride/

SEE ALSO

    dss, dssp
    '''
    from subprocess import Popen, PIPE
    import tempfile, os

    state, quiet = int(state), int(quiet)

    ss_map = {
        'C': 'L',
        'B': 'S',
        'b': 'S',
        'E': 'S',
        'T': 'L',
        'G': 'H',
        'H': 'H',
    }
    tmpfilepdb = tempfile.mktemp('.pdb')
    ss_dict = dict()
    for model in cmd.get_object_list(selection):
        cmd.save(tmpfilepdb, '%s and (%s)' % (model, selection), state)
        try:
            process = Popen([exe, tmpfilepdb], stdout=PIPE)
        except OSError:
            print 'Error: Cannot execute exe=' + exe
            raise CmdException
        for line in process.stdout:
            if not line.startswith('ASG'):
                continue
            chain = line[9].strip('-')
            resi = line[11:16].strip()
            ss = line[24]
            ss_dict[model,chain,resi] = ss
    os.remove(tmpfilepdb)
    _common_ss_alter(selection, ss_dict, ss_map, raw)

def set_phipsi(selection, phi=None, psi=None, state=1, quiet=1):
    '''
DESCRIPTION

    Set phi/psi angles for all residues in selection.

SEE ALSO

    phi_psi, cmd.get_phipsi, set_dihedral, DynoPlot
    '''
    for idx in cmd.index('byca (' + selection + ')'):
        x = cmd.index('((%s`%d) extend 2 and name C+N+CA)' % idx)
        if len(x) != 5 or x[2] != idx:
            print ' Warning: set_phipsi: missing atoms (%s`%d)' % idx
            continue
        if phi is not None:
            cmd.set_dihedral(x[0], x[1], x[2], x[3], phi, state, quiet)
        if psi is not None:
            cmd.set_dihedral(x[1], x[2], x[3], x[4], psi, state, quiet)

cmd.extend('split_chains', split_chains)
cmd.extend('rmsf2b', rmsf2b)
cmd.extend('set_sequence', set_sequence)
cmd.extend('alphatoall', alphatoall)
cmd.extend('mse2met', mse2met)
cmd.extend('polyala', polyala)
cmd.extend('stub2ala', stub2ala)
cmd.extend('remove_alt', remove_alt)
cmd.extend('dssp', dssp)
cmd.extend('stride', stride)
cmd.extend('set_phipsi', set_phipsi)

# tab-completion of arguments
cmd.auto_arg[0].update({
    'split_chains'   : cmd.auto_arg[0]['zoom'],
    'rmsf2b'         : cmd.auto_arg[0]['zoom'],
    'mse2met'        : cmd.auto_arg[0]['zoom'],
    'polyala'        : cmd.auto_arg[0]['zoom'],
    'stub2ala'       : cmd.auto_arg[0]['zoom'],
    'remove_alt'     : cmd.auto_arg[0]['zoom'],
    'dssp'           : cmd.auto_arg[0]['zoom'],
    'stride'         : cmd.auto_arg[0]['zoom'],
    'set_phipsi'     : cmd.auto_arg[0]['zoom'],
})
cmd.auto_arg[1].update({
    'set_sequence'   : cmd.auto_arg[0]['zoom'],
})

# vi: expandtab:smarttab
