'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

import sys

from pymol import cmd, CmdException

_auto_arg0_zoom = cmd.auto_arg[0]['zoom']


def _assert_package_import():
    if not __name__.endswith('.editing'):
        raise CmdException("Must do 'import psico.editing' instead of 'run ...'")


def split(operator, selection, prefix='entity', *, _self=cmd):
    '''
DESCRIPTION

    Create a single object for each entity in selection, defined by operator
    (e.g. bymolecule, bysegment, ...). Returns the number of created objects.
    '''
    _self.disable(' '.join(_self.get_object_list(selection)))
    tmp = _self.get_unused_name('_')
    _self.create(tmp, selection)

    r = 0
    while _self.count_atoms(tmp) > 0:
        name = _self.get_unused_name(prefix)
        _self.extract(name, operator + ' first model ' + tmp)
        r += 1

    _self.delete(tmp)
    return r


def split_molecules(selection='(all)', prefix='mol_', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Create a single object for each molecule (covalently connected entity) in
    selection (ignores solvent).

SEE ALSO

    split_chains, split_states
    '''
    quiet = int(quiet)
    r = split('bm.', '(%s) and not solvent' % selection, prefix, _self=_self)
    if not quiet:
        print(' Found %d non-solvent molecules' % r)


def split_chains(selection='(all)', prefix=None, *, _self=cmd):
    '''
DESCRIPTION

    Create a single object for each chain in selection

SEE ALSO

    split_states
    '''
    count = 0
    models = _self.get_object_list(selection)
    for model in models:
        for chain in _self.get_chains('(%s) and model %s' % (selection, model)):
            count += 1
            if not prefix:
                name = '%s_%s' % (model, chain)
            else:
                name = '%s%04d' % (prefix, count)
            _self.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
        _self.disable(model)


@cmd.extend
def split_segis(selection='all', prefix=None, *, _self=cmd):
    '''
DESCRIPTION

    Create a single object for each segi in selection

SEE ALSO

    split_chains
    '''
    from . import querying
    count = 0
    for model in _self.get_object_list(selection):
        sele = f"({selection}) & model {model}"
        for segi in set(querying.iterate_to_list(sele, "segi", _self=_self)):
            count += 1
            name = f'{prefix}{count:04d}' if prefix else f'{model}_{segi}'
            _self.create(name, f'({sele}) & segi "{segi}"')
        _self.disable(model)


def rmsf2b(selection='all', linearscale=1.0, var='b', quiet=1, *, _self=cmd):
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
    _assert_package_import()
    from . import querying
    from numpy import sqrt, pi
    linearscale = float(linearscale)
    coords = querying.get_ensemble_coords(selection, _self=_self)
    n_states, n_atoms, _ = coords.shape
    if n_atoms == 0 or n_states < 2:
        raise CmdException('not enough atoms or states')
    u_sq = coords.var(0).sum(1)  # var over states, sum over x,y,z
    b_array = sqrt(u_sq) * linearscale if linearscale > 0.0 \
            else 8 * pi**2 * u_sq
    _self.alter(selection, var + ' = next(b_iter)', space={'b_iter': iter(b_array), 'next': next})
    if not int(quiet):
        print(' Average RMSF: %.2f' % (sqrt(u_sq).mean()))
    return b_array


def set_sequence(sequence, selection='all', start=1, *, _self=cmd):
    '''
DESCRIPTION

    Alters the residue names according to given sequence

ARGUMENTS

    sequence = string: amino acid sequence in one-letter code

    selection = string: atom selection {default: all}

    start = int: residue number to start from {default: 1}
    '''
    import re
    _assert_package_import()
    from . import three_letter
    sequence = re.sub(r'\s+', '', sequence)
    start = int(start)
    for i, aa in enumerate(sequence):
        _self.alter('(%s) and resi %d' % (selection, i + start),
                'resn=' + repr(three_letter.get(aa.upper(), 'UNK')))


def alphatoall(selection='polymer', properties='b', operator='byca', quiet=1, *, _self=cmd):
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
    _self.iterate('%s (%s)' % (operator, selection), 'props[model,segi,chain,resi] = ' + properties,
            space=space)
    _self.alter(selection,
            properties + ' = props.get((model,segi,chain,resi), ' + properties + ')',
            space=space)
    if not int(quiet):
        print(' Modified %d residues' % (len(space['props'])))


def mse2met(selection='all', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Mutate selenomethionine to methionine
    '''
    quiet = int(quiet)
    x = _self.alter('(%s) and MSE/SE' % selection, 'name="SD";elem="S"')
    _self.flag('ignore', '(%s) and MSE/' % (selection), 'clear')
    _self.alter('(%s) and MSE/' % selection, 'resn="MET";type="ATOM"')
    if not quiet:
        print('Altered %d MSE residues to MET' % (x))
    _self.sort()


def polyala(selection='all', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Mutate any residue to Alanine (except Glycines)

SEE ALSO

    stub2ala
    '''
    quiet = int(quiet)
    _self.remove('polymer and (%s) and not name C+N+O+CA+CB+OXT' % (selection))
    _self.alter('polymer and (%s) and not resn GLY' % (selection), 'resn = "ALA"')
    _self.sort()


def stub2ala(selection='all', quiet=1, *, _self=cmd):
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
    _self.iterate('(%s) and polymer and (not hydro) and (not name C+N+O+OXT)' % (selection),
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
            _self.remove(key_str_cg)
            name_tuple = ('CA', 'CB')
        lookslike_resn = lookslike.get(name_tuple, resn)
        if lookslike_resn != resn:
            if not quiet:
                print('Altering %s to %s' % (key_str, lookslike_resn))
            _self.alter(key_str, 'resn = %s' % (repr(lookslike_resn)))
    _self.sort()


def remove_alt(selection='all', keep='first', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Remove alternative location atoms.

USAGE

    remove_alt [selection [, keep]]

ARGUMENTS

    selection = string: atom selection

    keep = string: AltLoc to keep, or 'first' to keep the first observed AltLoc {default: first}
    '''
    if keep == "first":
        return remove_alt_keep_first(selection, quiet=quiet, _self=_self)

    if len(keep) != 1:
        raise CmdException(
            f"keep must be 'first' or a single letter, got {keep!r}")

    _self.remove('(%s) and not alt +%s' % (selection, keep), quiet=int(quiet))
    _self.alter(selection, '(alt,q)=("",1.0)')
    _self.sort()


def remove_alt_keep_first(selection='*', *, quiet=1, _self=cmd):
    '''
    Remove alternative location atoms, keep the first observed.
    '''
    alts = {}
    expr = "(alt, q) = callback((model, segi, chain, resi, resn, name), alt)"

    def callback(namekey, alt):
        return ("#", 0.) if alts.setdefault(namekey, alt) != alt else ("", 1.)

    tmpsele = _self.get_unused_name("_sele")
    _self.select(tmpsele, selection)
    try:
        _self.alter(tmpsele, expr, space={"callback": callback})
        _self.remove(f"{tmpsele} & not alt ''", quiet=quiet)
    finally:
        _self.delete(tmpsele)


def _common_ss_alter(selection, ss_dict, ss_map, raw='', *, _self=cmd):
    '''
DESCRIPTION

    Shared code of 'dssp' and 'stride' functions.
    '''
    if raw != 'ss':
        _self.alter(selection, 'ss = ss_map.get(ss_dict.get((model,chain,resi)), "")',
                space={'ss_dict': ss_dict, 'ss_map': ss_map})
    if raw != '':
        _self.alter(selection, raw + ' = ss_dict.get((model,chain,resi), "")',
                space={'ss_dict': ss_dict})
    _self.cartoon('auto', selection)
    _self.rebuild(selection, 'cartoon')


def dssp(selection='(all)', exe='', raw='', state=-1, quiet=1, *, _self=cmd):
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
    import subprocess, re
    import tempfile, os

    state, quiet = int(state), int(quiet)

    if exe == '':
        _assert_package_import()
        from . import which
        exe = which('dsspcmbi', 'dssp', 'dssp-2', 'mkdssp')

    args = [exe]

    try:
        version_stdout = subprocess.check_output([exe, "--version"])
        version = int(re.findall(br"version (\d+)", version_stdout)[0])
    except Exception as ex:
        print(f" dssp-Warning: {ex}")
        version = 0

    if version >= 4:
        args += ["--output-format", "dssp"]

    ss_map = {
        'B': 'S',  # residue in isolated beta-bridge
        'E': 'S',  # extended strand, participates in beta ladder
        'T': 'L',  # hydrogen bonded turn
        'G': 'H',  # 3-helix (3/10 helix)
        'H': 'H',  # alpha helix
        'I': 'H',  # 5 helix (pi helix)
        'S': 'L',  # bend
        ' ': 'L',  # loop or irregular
    }
    tmpfilepdb = tempfile.mktemp('.pdb')
    ss_dict = dict()
    for model in _self.get_object_list(selection):
        cmd.multisave(tmpfilepdb,
                '%s and (%s)' % (model, selection), state, _self=_self)
        try:
            process = subprocess.run(args + [tmpfilepdb],
                    capture_output=True,
                    universal_newlines=True)
        except OSError as ex:
            raise CmdException('Cannot execute exe=' + exe) from ex
        if process.returncode != 0:
            raise CmdException('dssp failed: ' + process.stderr.strip())
        line_it = iter(process.stdout.splitlines())
        for line in line_it:
            if line.startswith('  #  RESIDUE'):
                break
        for line in line_it:
            resi = line[5:11].strip()
            chain = line[11].strip()
            ss = line[16]
            ss_dict[model, chain, resi] = ss
    os.remove(tmpfilepdb)
    _common_ss_alter(selection, ss_dict, ss_map, raw, _self=_self)


def stride(selection='(all)', exe='stride', raw='', state=-1, quiet=1, *, _self=cmd):
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
    for model in _self.get_object_list(selection):
        _self.save(tmpfilepdb, '%s and (%s)' % (model, selection), state)
        try:
            process = Popen([exe, tmpfilepdb], stdout=PIPE,
                    universal_newlines=True)
        except OSError as ex:
            raise CmdException('Cannot execute exe=' + exe) from ex
        for line in process.stdout:
            if not line.startswith('ASG'):
                continue
            chain = line[9].strip('-')
            resi = line[11:16].strip()
            ss = line[24]
            ss_dict[model, chain, resi] = ss
    os.remove(tmpfilepdb)
    _common_ss_alter(selection, ss_dict, ss_map, raw, _self=_self)


def dss_promotif(selection='all', exe='', raw='', state=-1, quiet=1, *, _self=cmd):
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
        _assert_package_import()
        from . import which
        motifdir = os.environ.get('motifdir')
        exe = which('p_sstruc3', 'p_sstruc2', 'promotif.scr',
                path=[motifdir] if motifdir else None)

    tmpdir = tempfile.mkdtemp()
    tmpfilepdb = os.path.join(tmpdir, 'xxxx.pdb')
    tmpfilesst = os.path.join(tmpdir, 'xxxx.sst')
    ss_dict = dict()

    try:
        for model in _self.get_object_list(selection):
            _self.save(tmpfilepdb, 'model %s and (%s)' % (model, selection), state)

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
                    ss_dict[model, chain, resi] = ss

            os.remove(tmpfilesst)
    except OSError as ex:
        raise CmdException('Cannot execute exe=' + exe) from ex
    finally:
        shutil.rmtree(tmpdir)
    _common_ss_alter(selection, ss_dict, ss_map, raw, _self=_self)


def sst(selection='(all)', raw='', state=-1, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Secondary structure assignment with SST.
    https://lcb.infotech.monash.edu/sst/sst_submission_form.html

SEE ALSO

    dss, dssp, stride
    '''
    import ssl
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

    for model in _self.get_object_list(selection):
        pdbstr = _self.get_pdbstr('%s & guide & (%s)' % (model, selection), state)

        body = '\r\n'.join([
            '--' + boundary,
            'Content-Disposition: form-data; name="sstverstr"',
            'Content-Type: text/plain',
            '',
            'v1.6',
            '--' + boundary,
            'Content-Disposition: file; name="pdb_file"; filename="abc.pdb"',
            'Content-Type: text/plain',
            '',
            pdbstr,
            '--' + boundary + '--',
            '',
        ])

        body = body.encode('ascii', 'ignore')

        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE

        try:
            request = urllib2.Request(
                data=body,
                url='https://lcb.infotech.monash.edu/sst/sst_formaction.php')
            request.add_header('User-agent', 'PyMOL ' + cmd.get_version()[0] + ' ' +
                    sys.platform)
            request.add_header('Content-type', 'multipart/form-data; boundary=%s' % boundary)
            request.add_header('Content-length', len(body))
            lines = urllib2.urlopen(request, context=ctx).readlines()
        except urllib2.URLError as ex:
            raise CmdException('URL request failed') from ex

        # Example output
        """
....................RESIDUE-LEVEL ASSIGNMENT DETAILS......................
.
! Chain Residue     SSType
- A     M (MET1)    C
- A     G (GLY76)   C
.
..........................................................................
...................................LEGEND.................................
        """

        lines = (line.decode('ascii', 'ignore') for line in lines)
        lines = iter(lines)

        for line in lines:
            if line.startswith('....................RESIDUE-LEVEL'):
                break
        else:
            if not quiet:
                print(' Warning: SST assignment failed')
            return

        next(lines)

        for line in lines:
            if 'LEGEND........' in line:
                break
            if not line.startswith("- "):
                continue
            chain = line[2:6].rstrip()
            # resn = line[11:14].strip()
            resi = line[14:19].rstrip(") ")
            ss = line[20]
            ss_dict[model, chain, resi] = ss

    _common_ss_alter(selection, ss_dict, ss_map, raw, _self=_self)


def set_phipsi(selection, phi=None, psi=None, state=1, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Set phi/psi angles for all residues in selection.

SEE ALSO

    phi_psi, cmd.get_phipsi, set_dihedral, DynoPlot
    '''
    for idx in _self.index('byca (' + selection + ')'):
        x = _self.index('((%s`%d) extend 2 and name C+N+CA)' % idx)
        if len(x) != 5 or x[2] != idx:
            print(' Warning: set_phipsi: missing atoms (%s`%d)' % idx)
            continue
        if phi is not None:
            _self.set_dihedral(x[0], x[1], x[2], x[3], phi, state, quiet)
        if psi is not None:
            _self.set_dihedral(x[1], x[2], x[3], x[4], psi, state, quiet)


def update_identifiers(target, source, identifiers='segi chain resi',
        match='align', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Transfers segi, chain, and resi identifiers from one selection to another.
    This works by mapping old to new identifiers and alters also not aligned
    atoms (works if any other atom from the same residue got aligned).
    '''
    from .fitting import MatchMaker

    with MatchMaker(target, source, match, _self=_self) as mm:
        key = '(' + ','.join(identifiers.split()) + ',)'
        tkeys, skeys = [], []
        _self.iterate(mm.mobile, 'tkeys.append(%s)' % (key), space=locals())
        _self.iterate(mm.target, 'skeys.append(%s)' % (key), space=locals())
        t2s = dict(zip(tkeys, skeys))
        _self.alter(target, '%s = t2s.get(%s, %s)' % (key, key, key), space=locals())


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
    'split_chains': _auto_arg0_zoom,
    'split_molecules': _auto_arg0_zoom,
    'split_segis': _auto_arg0_zoom,
    'rmsf2b': _auto_arg0_zoom,
    'mse2met': _auto_arg0_zoom,
    'polyala': _auto_arg0_zoom,
    'stub2ala': _auto_arg0_zoom,
    'remove_alt': _auto_arg0_zoom,
    'dssp': _auto_arg0_zoom,
    'stride': _auto_arg0_zoom,
    'dss_promotif': _auto_arg0_zoom,
    'sst': _auto_arg0_zoom,
    'set_phipsi': _auto_arg0_zoom,
})
cmd.auto_arg[1].update({
    'set_sequence': _auto_arg0_zoom,
})

# vi: expandtab:smarttab
