'''
(c) 2010-2012 Thomas Holder

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException

def join_states(name, selection='all', discrete=-1, zoom=0, quiet=1):
    '''
DESCRIPTION

    The reverse of split_states

ARGUMENTS

    name = string: name of object to create or modify
 
    selection = string: atoms to include in the new object

    discrete = -2: match atoms by sequence alignment
    discrete = -1: Assume identical input objects (matching all atom
        identifiers) but also check for missing atoms and only include atoms
        that are present in all input objects {default}
    discrete = 0: Assume identical input objects
    discrete = 1: Input object may be (totally) different
    '''
    discrete, quiet = int(discrete), int(quiet)
    if discrete == -2:
        from .selecting import wait_for
        aln_obj = cmd.get_unused_name('_')
    models = cmd.get_object_list('(' + selection + ')')
    for i in range(len(models)):
        if discrete == -1 and i > 0:
            cmd.remove('(%s) and not (alt A+ and (%s) in (%s))' % (name, name, models[i]))
            cmd.create(name, '(%s) in (%s)' % (models[i], name), 1, i+1, 0, 0, quiet)
        elif discrete == -2 and i > 0:
            cmd.align(models[i], name, cycles=0, transform=0, object=aln_obj)
            wait_for(aln_obj)
            cmd.remove('(%s) and not (%s)' % (name, aln_obj))
            cmd.create(name, name, 1, i+1, 0, 0, quiet)
            cmd.update(name, '(%s) and (%s)' % (models[i], aln_obj), i+1, 1, 0, quiet)
            cmd.delete(aln_obj)
        else:
            cmd.create(name, models[i], 1, i+1, discrete == 1, 0, quiet)
    if int(zoom):
        cmd.zoom(name, state=0)

sidechaincenteratoms = {
    'GLY': ('CA',),
    'ALA': ('CB',),
    'VAL': ('CG1', 'CG2'),
    'ILE': ('CD1',),
    'LEU': ('CD1', 'CD2'),
    'SER': ('OG',),
    'THR': ('OG1', 'CG2'),
    'ASP': ('OD1', 'OD2'),
    'ASN': ('OD1', 'ND2'),
    'GLU': ('OE1', 'OE2'),
    'GLN': ('OE1', 'NE2'),
    'LYS': ('NZ',),
    'ARG': ('NE', 'NH1', 'NH2'),
    'CYS': ('SG',),
    'MET': ('SD',),
    'MSE': ('SE',),
    'PHE': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
    'TYR': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'),
    'TRP': ('CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3'),
    'HIS': ('CG', 'ND1', 'CD2', 'CE1', 'NE2'),
    'PRO': ('CB', 'CG', 'CD'),
}

sidechaincentermethods = ['bahar1996', 'centroid']

def sidechaincenters(object='scc', selection='all', method='bahar1996', name='PS1'):
    '''
DESCRIPTION

    Creates an object with sidechain representing pseudoatoms for each residue
    in selection.

    Two methods are available:
    (1) Sidechain interaction centers as defined by Bahar and Jernigan 1996
        http://www.ncbi.nlm.nih.gov/pubmed/9080182
    (2) Sidechain centroids, the pseudoatom is the centroid of all atoms except
        hydrogens and backbone atoms (N, C and O).

NOTE

    With method "bahar1996", if a residue has all relevant sidechain center
    atoms missing (for example a MET without SD), it will be missing in the
    created pseudoatom object.

    With method "centroid", if you want to exclude C-alpha atoms from
    sidechains, modify the selection like in this example:

    sidechaincenters newobject, all and (not name CA or resn GLY), method=2

USAGE

    sidechaincenters object [, selection [, method ]]

ARGUMENTS

    object = string: name of object to create

    selection = string: atoms to consider {default: (all)}

    method = string: bahar1996 or centroid {default: bahar1996}

    name = string: atom name of pseudoatoms {default: PS1}

SEE ALSO

    pseudoatom
    '''
    from chempy import Atom, cpv, models

    atmap = dict()
    if method in ['bahar1996', '1', 1]:
        modelAll = cmd.get_model('(%s) and resn %s' % (selection, '+'.join(sidechaincenteratoms)))
        for at in modelAll.atom:
            if at.name in sidechaincenteratoms[at.resn]:
                atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)
    elif method in ['centroid', '2', 2]:
        modelAll = cmd.get_model('(%s) and polymer and not (hydro or name C+N+O)' % selection)
        for at in modelAll.atom:
            atmap.setdefault((at.segi, at.chain, at.resn, at.resi), []).append(at)
    else:
        print('Error: unknown method:', method)
        raise CmdException

    model = models.Indexed()
    for centeratoms in atmap.values():
        center = cpv.get_null()
        for at in centeratoms:
            center = cpv.add(center, at.coord)
        center = cpv.scale(center, 1./len(centeratoms))
        atom = Atom()
        atom.coord = center
        atom.index = model.nAtom + 1
        atom.name = name
        for key in ['resn','chain','resi','resi_number','hetatm','ss','segi']:
            atom.__dict__[key] = at.__dict__[key]
        model.add_atom(atom)
    model.update_index()
    if object in cmd.get_object_list():
        cmd.delete(object)
    cmd.load_model(model, object)
    return model

def ramp_levels(name, levels, quiet=1):
    '''
DESCRIPTION

    Changes the slot levels of a ramp.

SEE ALSO

    ramp_new, isolevel
    '''
    quiet = int(quiet)
    if cmd.is_string(levels):
        levels = cmd.safe_list_eval(levels)

    try:
        position = cmd.get_names('all').index(name)

        odata = cmd.get_session(name, 1, 1, 0, 0)['names'][0]
        if odata[4] != 8:
            raise TypeError('not a ramp')
        data = odata[5]
    except:
        print(' Error: Get session data for ramp "%s" failed' % (name))
        raise CmdException

    if len(levels) != len(data[3]):
        print(' Error: number of levels must agree with existing object')
        raise CmdException

    map_name = data[6]
    colors = [data[4][i:i+3] for i in range(0, len(data[4]), 3)]
    cmd.ramp_new(name, map_name, levels, colors, quiet=quiet)

    # restore original position
    if position == 0:
        cmd.order(name, location='top')
    else:
        cmd.order(cmd.get_names('all')[position-1] + ' ' + name)

def pdb2pqr(name, selection='all', ff='amber', debump=1, opt=1, assignonly=0,
        ffout='', ph=None, neutraln=0, neutralc=0, state=-1, preserve=0,
        exe='pdb2pqr', quiet=1):
    '''
DESCRIPTION

    Creates a new molecule object from a selection and adds missing atoms,
    assignes charges and radii using PDB2PQR.

    http://www.poissonboltzmann.org/pdb2pqr/

USAGE

    pdb2pqr name [, selection [, ff [, debump [, opt [, assignonly [, ffout [,
        ph [, neutraln [, neutralc [, state [, preserve ]]]]]]]]]]]

ARGUMENTS

    name = string: name of object to create or modify

    selection = string: atoms to include in the new object {default: all}

    ff = string: forcefield {default: amber}
    '''
    import os, tempfile, subprocess, shutil

    debump, opt, assignonly = int(debump), int(opt), int(assignonly)
    neutraln, neutralc = int(neutraln), int(neutralc)
    state, preserve, quiet = int(state), int(preserve), int(quiet)

    args = [cmd.exp_path(exe), '--ff=' + ff, '--chain']
    if not debump:
        args.append('--nodebump')
    if not opt:
        args.append('--noopt')
    if assignonly:
        args.append('--assign-only')
    if ffout:
        args.append('--ffout=' + ffout)
    if ph is not None:
        args.append('--with-ph=%f' % float(ph))
    if neutraln:
        args.append('--neutraln')
    if neutralc:
        args.append('--neutralc')
    if not quiet:
        args.append('--verbose')

    tmpdir = tempfile.mkdtemp()
    infile = os.path.join(tmpdir, 'in.pdb')
    outfile = os.path.join(tmpdir, 'out.pqr')
    args.extend([infile, outfile])

    try:
        cmd.save(infile, selection, state)

        p = subprocess.Popen(args, cwd=tmpdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)

        print(p.communicate()[0].rstrip())

        if p.returncode != 0:
            raise CmdException('%s failed with exit status %d' % (args[0], p.returncode))

        remark5 = [L for L in open(outfile) if L.startswith('REMARK   5')]
        if remark5:
            print(''.join(remark5))

        cmd.load(outfile, name)
    except OSError:
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if not preserve:
            shutil.rmtree(tmpdir)
        elif not quiet:
            print(' Notice: not deleting', tmpdir)

    if not quiet:
        print(' pdb2pqr: done')

def corina(name, selection, exe='corina', state=-1, preserve=0, quiet=1):
    '''
DESCRIPTION

    Run corina on selection and load as new object
    '''
    import os, tempfile, subprocess, shutil
    from . import querying

    state, preserve, quiet = int(state), int(preserve), int(quiet)

    if state < 1:
        state = querying.get_selection_state(selection)

    tmpdir = tempfile.mkdtemp()
    infile = os.path.join(tmpdir, 'in.sdf')
    outfile = os.path.join(tmpdir, 'out.sdf')
    trcfile = os.path.join(tmpdir, 'corina.trc')

    args = [exe, infile, outfile]
    stderr = ''

    try:
        cmd.save(infile, selection, state)
        stdout, stderr = subprocess.Popen(args, cwd=tmpdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE).communicate()

        if not quiet and stdout.strip():
            print(stdout)

        trclines = open(trcfile).readlines()
        trcerror = [line for line in trclines if 'ERROR' in line]
        if trcerror:
            raise CmdException("corina failed: " + trcerror[0].strip())

        if not os.path.exists(outfile):
            raise CmdException("corina failed: " + stderr.strip())

        cmd.load(outfile, name)
    except OSError:
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if not preserve:
            shutil.rmtree(tmpdir)
        elif not quiet:
            print(' Notice: not deleting', tmpdir)

    if not quiet:
        print(' corina: done')

def prepwizard(name, selection='all', options='', state=-1,
        preserve=0, exe='$SCHRODINGER/utilities/prepwizard', quiet=1):
    '''
DESCRIPTION

    Run the SCHRODINGER Protein Preparation Wizard. Builds missing side
    chains and converts MSE to MET. Other non-default options need to be
    passed with the "options=" argument.

USAGE

    prepwizard name [, selection [, options [, state ]]]

ARGUMENTS

    name = str: name of object to create

    selection = str: atoms to send to the wizard {default: all}

    options = str: additional command line options for prepwizard

    state = int: object state {default: -1 (current)}
    '''
    import os, tempfile, subprocess, shutil, shlex

    state, preserve, quiet = int(state), int(preserve), int(quiet)

    exe = cmd.exp_path(exe)
    if not os.path.exists(exe):
        if 'SCHRODINGER' not in os.environ:
            print(' Warning: SCHRODINGER environment variable not set')
        raise CmdException('no such script: ' + exe)

    args = [exe, '-mse', '-fillsidechains', '-WAIT']

    if options:
        if cmd.is_string(options):
            options = shlex.split(options)
        args.extend(options)

    tmpdir = tempfile.mkdtemp()
    infile = 'in.pdb'
    outfile = 'out.mae'
    args.extend([infile, outfile])

    try:
        cmd.save(os.path.join(tmpdir, infile), selection, state)

        p = subprocess.Popen(args, cwd=tmpdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)

        print(p.communicate()[0].rstrip())

        if p.returncode != 0:
            logfile = os.path.join(tmpdir, 'in.log')
            if os.path.exists(logfile):
                with open(logfile) as handle:
                    print(handle.read())

            raise CmdException('%s failed with exit status %d' % (args[0], p.returncode))

        cmd.load(os.path.join(tmpdir, outfile), name)
    except OSError:
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if not preserve:
            shutil.rmtree(tmpdir)
        elif not quiet:
            print(' Notice: not deleting', tmpdir)

    if not quiet:
        print(' prepwizard: done')

def fiber(seq, num=4, name='', rna=0, single=0, repeats=0,
        preserve=0, exe='$X3DNA/bin/fiber', quiet=1):
    '''
DESCRIPTION

    Run X3DNA's "fiber" tool.

    For the list of structure identification numbers, see for example:
    http://xiang-jun.blogspot.com/2009/10/fiber-models-in-3dna.html

USAGE

    fiber seq [, num [, name [, rna [, single ]]]]

ARGUMENTS

    seq = str: single letter code sequence or number of repeats for
    repeat models.

    num = int: structure identification number {default: 4}

    name = str: name of object to create {default: random unused name}

    rna = 0/1: 0=DNA, 1=RNA {default: 0}

    single = 0/1: 0=double stranded, 1=single stranded {default: 0}

EXAMPLES

    # environment (this could go into ~/.pymolrc or ~/.bashrc)
    os.environ["X3DNA"] = "/opt/x3dna-v2.3"

    # B or A DNA from sequence
    fiber CTAGCG
    fiber CTAGCG, 1, ADNA

    # double or single stranded RNA from sequence
    fiber AAAGGU, name=dsRNA, rna=1
    fiber AAAGGU, name=ssRNA, rna=1, single=1

    # poly-GC Z-DNA repeat model with 10 repeats
    fiber 10, 15
    '''
    import os, tempfile, subprocess, shutil

    rna, single = int(rna), int(single)
    preserve, quiet = int(preserve), int(quiet)

    if 'X3DNA' not in os.environ:
        raise CmdException('please set X3DNA environment variable')

    args = [cmd.exp_path(exe), '-seq=' + seq, '-' + str(num)]

    if rna:
        args.append('-rna')
    if single:
        args.append('-single')

    if not name:
        name = cmd.get_unused_name('fiber-' + str(num) + '_')

    tmpdir = tempfile.mkdtemp()
    outfile = os.path.join(tmpdir, 'out.pdb')
    args.append(outfile)

    if seq.endswith('help'):
        # call fiber with no arguments to get help page
        args[1:] = []

    try:
        p = subprocess.Popen(args, cwd=tmpdir,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)

        stdout, _ = p.communicate(str(int(repeats) or seq))

        if not quiet:
            print(stdout)

        if len(args) != 1:
            if p.returncode != 0:
                raise CmdException('Returned non-zero status: ' + str(args))

            cmd.load(outfile, name, quiet=quiet)

    except OSError:
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if not preserve:
            shutil.rmtree(tmpdir)
        elif not quiet:
            print(' Notice: not deleting', tmpdir)

if 'join_states' not in cmd.keyword:
    cmd.extend('join_states', join_states)
cmd.extend('sidechaincenters', sidechaincenters)
cmd.extend('ramp_levels', ramp_levels)
cmd.extend('pdb2pqr', pdb2pqr)
cmd.extend('corina', corina)
cmd.extend('prepwizard', prepwizard)
cmd.extend('fiber', fiber)

cmd.auto_arg[0].update([
    ('ramp_levels', [lambda: cmd.Shortcut(cmd.get_names_of_type('object:')), 'ramp object', '']),
])
cmd.auto_arg[1].update([
    ('join_states'          , cmd.auto_arg[0]['zoom']),
    ('sidechaincenters'     , cmd.auto_arg[0]['zoom']),
    ('pdb2pqr'              , cmd.auto_arg[0]['zoom']),
    ('corina'               , cmd.auto_arg[0]['zoom']),
    ('prepwizard'           , cmd.auto_arg[0]['zoom']),
])
cmd.auto_arg[2].update([
    ('sidechaincenters', [cmd.Shortcut(sidechaincentermethods), 'method', '']),
    ('pdb2pqr', [cmd.Shortcut(['amber', 'charmm', 'parse', 'tyl06']), 'forcefield', '']),
])

# vi: ts=4:sw=4:smarttab:expandtab
