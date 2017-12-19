'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import absolute_import
from __future__ import with_statement

import sys

from pymol import cmd, CmdException

def split(operator, selection, prefix='entity'):
    '''
DESCRIPTION

    Create a single object for each entity in selection, defined by operator
    (e.g. bymolecule, bysegment, ...). Returns the number of created objects.
    '''
    cmd.disable(' '.join(cmd.get_object_list('(' + selection + ')')))
    tmp = cmd.get_unused_name('_')
    cmd.create(tmp, selection)

    r = 0
    while cmd.count_atoms(tmp) > 0:
        name = cmd.get_unused_name(prefix)
        cmd.extract(name, operator + ' first model ' + tmp)
        r += 1

    cmd.delete(tmp)
    return r

def split_molecules(selection='(all)', prefix='mol_', quiet=1):
    '''
DESCRIPTION

    Create a single object for each molecule (covalently connected entity) in
    selection (ignores solvent).

SEE ALSO

    split_chains, split_states
    '''
    quiet = int(quiet)
    r = split('bm.', '(%s) and not solvent' % selection, prefix)
    if not quiet:
        print(' Found %d non-solvent molecules' % r)

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

def rmsf2b(selection='all', linearscale=1.0, var='b', quiet=1):
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
    from numpy import asfarray, sqrt, pi
    linearscale = float(linearscale)
    n_atoms = cmd.count_atoms(selection)
    n_states = cmd.count_states(selection)
    if n_atoms == 0 or n_states < 2:
        print(' Error: not enough atoms or states')
        raise CmdException
    coords = []
    cmd.iterate_state(0, selection, 'coords.append((x,y,z))', atomic=0,
            space={'coords': coords})
    coords = asfarray(coords).reshape((cmd.count_states(selection), -1, 3))
    u_sq = coords.var(0).sum(1) # var over states, sum over x,y,z
    b_array = sqrt(u_sq) * linearscale if linearscale > 0.0 \
            else 8 * pi**2 * u_sq
    cmd.alter(selection, var + ' = next(b_iter)', space={'b_iter': iter(b_array), 'next': next})
    if not int(quiet):
        print(' Average RMSF: %.2f' % (sqrt(u_sq).mean()))
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
        print(' Modified %d residues' % (len(space['props'])))

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
        print('Altered %d MSE residues to MET' % (x))
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
                print('Removing ' + str(key_str_cg))
            cmd.remove(key_str_cg)
            name_tuple = ('CA', 'CB')
        lookslike_resn = lookslike.get(name_tuple, resn)
        if lookslike_resn != resn:
            if not quiet:
                print('Altering %s to %s' % (key_str, lookslike_resn))
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
    cmd.cartoon('auto', selection)
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
    from .exporting import save_pdb_without_ter

    state, quiet = int(state), int(quiet)

    if exe == '':
        from . import which
        exe = which('dsspcmbi', 'dssp', 'dssp-2', 'mkdssp')
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
        save_pdb_without_ter(tmpfilepdb,
                '%s and (%s)' % (model, selection), state)
        try:
            process = Popen([exe, tmpfilepdb], stdout=PIPE)
        except OSError:
            print('Error: Cannot execute exe=' + exe)
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
            print('Error: Cannot execute exe=' + exe)
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

def dss_promotif(selection='all', exe='', raw='', state=-1, quiet=1):
    '''
DESCRIPTION

    Secondary structure assignment with PROMOTIF.
    http://www.rubic.rdg.ac.uk/~gail/#Software

SEE ALSO

    dss, dssp, stride
    '''
    from subprocess import Popen, PIPE
    import tempfile, os, shutil

    state, quiet = int(state), int(quiet)

    ss_map = {
        'B': 'S',
        'E': 'S',
        'H': 'H',
        'G': 'H',
    }

    exe = cmd.exp_path(exe)
    if not exe:
        from . import which
        motifdir = os.environ.get('motifdir')
        exe = which('p_sstruc3', 'p_sstruc2', 'promotif.scr',
                path=[motifdir] if motifdir else None)

    tmpdir = tempfile.mkdtemp()
    tmpfilepdb = os.path.join(tmpdir, 'xxxx.pdb')
    tmpfilesst = os.path.join(tmpdir, 'xxxx.sst')
    ss_dict = dict()

    try:
        for model in cmd.get_object_list('(' + selection + ')'):
            cmd.save(tmpfilepdb, 'model %s and (%s)' % (model, selection), state)

            process = Popen([exe, tmpfilepdb], cwd=tmpdir, stdin=PIPE)
            process.communicate(tmpfilepdb + os.linesep)

            with open(tmpfilesst) as handle:
                for line in handle:
                    if line.startswith(' num  seq.no'):
                        break
                for line in handle:
                    if not line.strip():
                        break
                    chain = line[6].strip('-')
                    resi = line[7:12].strip()
                    ss = line[23]
                    ss_dict[model,chain,resi] = ss

            os.remove(tmpfilesst)
    except OSError:
        print(' Error: Cannot execute exe=' + exe)
        raise CmdException
    finally:
        shutil.rmtree(tmpdir)
    _common_ss_alter(selection, ss_dict, ss_map, raw)

def sst(selection='(all)', raw='', state=-1, quiet=1):
    '''
DESCRIPTION

    Secondary structure assignment with SST.
    http://lcb.infotech.monash.edu.au/sstweb/

SEE ALSO

    dss, dssp, stride
    '''
    try:
        import urllib2
    except ImportError:
        import urllib.request as urllib2

    state, quiet = int(state), int(quiet)

    ss_map = {
        'C': 'L',
        'E': 'S',
        'G': 'H',
        'H': 'H',
        'I': 'H',
        'g': 'H',
        'h': 'H',
        'i': 'H',
        '3': 'L',
        '4': 'L',
        '5': 'L',
        'T': 'L',
        '-': 'L',
        '|': 'L',
        ':': 'H',
    }
    ss_dict = {}
    boundary = '192.168.1.80.500.9981.1375391031.267.10'

    for model in cmd.get_object_list(selection):
        pdbstr = cmd.get_pdbstr('%s & guide & (%s)' % (model, selection), state)

        body = '\r\n'.join([
            '--' + boundary,
            'Content-Disposition: file; name="pdb_file"; filename="abc.pdb"',
            'Content-Type: text/plain',
            '',
            pdbstr,
            '--' + boundary + '--',
            '',
        ])

        body = body.encode('ascii', 'ignore')

        try:
            request = urllib2.Request(
                    data=body, url=
                    'http://lcb.infotech.monash.edu.au/sstweb/formaction_pdbfile.php')
            request.add_header('User-agent', 'PyMOL ' + cmd.get_version()[0] + ' ' +
                    cmd.sys.platform)
            request.add_header('Content-type', 'multipart/form-data; boundary=%s' % boundary)
            request.add_header('Content-length', len(body))
            lines = urllib2.urlopen(request).readlines()
        except urllib2.URLError:
            print(' Error: URL request failed')
            raise CmdException

        if sys.version_info[0] > 2:
            lines = (line.decode('ascii', 'ignore') for line in lines)

        lines = iter(lines)

        for line in lines:
            if line.startswith('..........RESIDUE LEVEL'):
                break
        else:
            if not quiet:
                print(' Warning: SST assignment failed')
            return

        next(lines)

        for line in lines:
            if line.startswith('...................................END'):
                break
            chain = line[2].strip()
            resi = line[3:9].strip()
            ss = line[21]
            ss_dict[model,chain,resi] = ss

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
            print(' Warning: set_phipsi: missing atoms (%s`%d)' % idx)
            continue
        if phi is not None:
            cmd.set_dihedral(x[0], x[1], x[2], x[3], phi, state, quiet)
        if psi is not None:
            cmd.set_dihedral(x[1], x[2], x[3], x[4], psi, state, quiet)

def update_identifiers(target, source, identifiers='segi chain resi',
        match='align', quiet=1):
    '''
DESCRIPTION

    Transfers segi, chain, and resi identifiers from one selection to another.
    This works by mapping old to new identifiers and alters also not aligned
    atoms (works if any other atom from the same residue got aligned).
    '''
    from .fitting import matchmaker

    tmatched, smatched, tmp_names = matchmaker(target, source, match)

    key = '(' + ','.join(identifiers.split()) + ',)'
    tkeys, skeys = [], []
    cmd.iterate(tmatched, 'tkeys.append(%s)' % (key), space=locals())
    cmd.iterate(smatched, 'skeys.append(%s)' % (key), space=locals())
    t2s = dict(zip(tkeys, skeys))
    cmd.alter(target, '%s = t2s.get(%s, %s)' % (key, key, key), space=locals())

    for name in tmp_names:
        cmd.delete(name)

if 'split_chains' not in cmd.keyword:
    cmd.extend('split_chains', split_chains)
cmd.extend('split_molecules', split_molecules)
cmd.extend('rmsf2b', rmsf2b)
cmd.extend('set_sequence', set_sequence)
if 'alphatoall' not in cmd.keyword:
    cmd.extend('alphatoall', alphatoall)
if 'mse2met' not in cmd.keyword:
    cmd.extend('mse2met', mse2met)
cmd.extend('polyala', polyala)
cmd.extend('stub2ala', stub2ala)
cmd.extend('remove_alt', remove_alt)
cmd.extend('dssp', dssp)
cmd.extend('stride', stride)
cmd.extend('dss_promotif', dss_promotif)
cmd.extend('sst', sst)
cmd.extend('set_phipsi', set_phipsi)
cmd.extend('update_identifiers', update_identifiers)

# tab-completion of arguments
cmd.auto_arg[0].update({
    'split_chains'   : cmd.auto_arg[0]['zoom'],
    'split_molecules': cmd.auto_arg[0]['zoom'],
    'rmsf2b'         : cmd.auto_arg[0]['zoom'],
    'mse2met'        : cmd.auto_arg[0]['zoom'],
    'polyala'        : cmd.auto_arg[0]['zoom'],
    'stub2ala'       : cmd.auto_arg[0]['zoom'],
    'remove_alt'     : cmd.auto_arg[0]['zoom'],
    'dssp'           : cmd.auto_arg[0]['zoom'],
    'stride'         : cmd.auto_arg[0]['zoom'],
    'dss_promotif'   : cmd.auto_arg[0]['zoom'],
    'sst'            : cmd.auto_arg[0]['zoom'],
    'set_phipsi'     : cmd.auto_arg[0]['zoom'],
})
cmd.auto_arg[1].update({
    'set_sequence'   : cmd.auto_arg[0]['zoom'],
})

# vi: expandtab:smarttab
