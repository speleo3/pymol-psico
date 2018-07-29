'''
Normal mode calculation using external apps or libraries.

(c) 2011-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import print_function

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from pymol import cmd, stored, CmdException
from pymol import selector

def normalmodes_pdbmat(selection, cutoff=10.0, force=1.0, mass='COOR',
        first=7, last=10, choose='LOWE', substruct='RESI', blocksize=4,
        exe='pdbmat', diag_exe='diagrtb',
        prefix='mode', states=7, factor=-1, clean=1, quiet=1, async_=-1, **kwargs):
    '''
DESCRIPTION

    PDBMAT and DIAGRTB wrapper

    Runs "pdbmat" and "diagrtb" and generates perturbed models for modes
    "first" to "last". WARNING: May run for a long time.

    PDBMAT computes the mass-weighted second derivatives energy matrix,
    using Tirion's model, that is, an elastic network model (ENM).
    In such models, close particles (atoms) are linked by springs.

    http://ecole.modelisation.free.fr/modes.html

NOTES

    Only considers ATOM records, if your model contains MSE residues or
    ligands that you want to consider, prepare it like this:
    mse2met (all)                # for MSE residues
    alter (hetatm), type='ATOM'  # for any hetatm

ARGUMENTS

    selection = string: atom selection

    cutoff = float: interaction distance cutoff in angstroem {default: 10}

    force = float: interaction force constant {default: 1.0}

    mass = string: origin of mass values {default: COOR}

    first = int: first mode to create perturbed model {default: 7}

    last = int: last mode to create perturbed model {default: 10}

    choose = string: eigenvalues chosen {default: LOWE}

    substruct = string: type of substructuring {default: RESI}

    blocksize = int: nb of residues per block {default: 4}
    '''
    args = [selection, cutoff, force, mass,
            first, last, choose, substruct, blocksize,
            exe, diag_exe, prefix, states, factor, clean, quiet]

    quiet, async_ = int(quiet), int(kwargs.pop('async', async_))
    if async_ < 0:
        async_ = not quiet

    if not async_:
        return _normalmodes(*args)

    from pymol.wizard.message import Message
    wiz = Message(['normalmodes: please wait ...', ''], dismiss=0)
    cmd.set_wizard(wiz)

    import threading
    t = threading.Thread(target=_normalmodes, args=args + [wiz])
    t.setDaemon(1)
    t.start()

def _normalmodes(selection, cutoff, force, mass,
        first, last, choose, substruct, blocksize,
        exe, diag_exe, prefix, states, factor, clean, quiet, wiz=None):
    import tempfile, subprocess, os, shutil, sys
    from chempy import cpv

    cutoff, force = float(cutoff), float(force)
    first, last, blocksize = int(first), int(last), int(blocksize)
    clean, quiet = int(clean), int(quiet)

    exe = cmd.exp_path(exe)
    diag_exe = cmd.exp_path(diag_exe)
    tempdir = tempfile.mkdtemp()

    if not quiet:
        print(' normalmodes: Temporary directory is', tempdir)

    try:
        sele_name = cmd.get_unused_name('__pdbmat')
    except AttributeError:
        sele_name = '__pdbmat'

    try:
        filename = os.path.join(tempdir, 'mobile.pdb')
        commandfile = os.path.join(tempdir, 'pdbmat.dat')

        cmd.select(sele_name, '(%s) and not hetatm' % (selection))
        cmd.save(filename, sele_name)
        
        f = open(commandfile, 'w')
        f.write('''! pdbmat file
 Coordinate FILENAME        = %s
 MATRIx FILENAME            = matrix.sdijb
 INTERACtion DISTance CUTOF = %.3f
 INTERACtion FORCE CONStant = %.3f
 Origin of MASS values      = %s
 Output PRINTing level      =          0
 MATRix FORMat              =       BINA
''' % (filename, cutoff, force, mass.upper()))
        f.close()

        if not quiet:
            print(' normalmodes: running', exe, '...')
        if wiz is not None:
            wiz.message[1] = 'running pdbmat'
            cmd.refresh_wizard()
        process = subprocess.Popen([exe], cwd=tempdir,
                stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

        for line in process.stdout:
            if not quiet:
                sys.stdout.write(line)

        try:
            natoms = len(open(os.path.join(tempdir, 'pdbmat.xyzm')).readlines())
        except Exception as e:
            print(e)
            natoms = cmd.count_atoms(sele_name)

        if natoms != cmd.count_atoms(sele_name):
            print('Error: pdbmat did not recognize all atoms')
            raise CmdException

        commandfile = os.path.join(tempdir, 'diagrtb.dat')
        f = open(commandfile, 'w')
        f.write('''! diagrtb file
 MATRIx FILENAME            = matrix.sdijb
 COORdinates filename       = %s
 Eigenvector OUTPut filename= diagrtb.eigenfacs
 Nb of VECTors required     = %d
 EigeNVALues chosen         = %s
 Type of SUBStructuring     = %s
 Nb of residues per BLOck   = %d
 Origin of MASS values      = %s
 Temporary files cleaning   =       ALL
 MATRix FORMat              =       BINA
 Output PRINting level      =          0
''' % (filename, last, choose.upper(), substruct.upper(), blocksize, mass.upper()))
        f.close()

        exe = diag_exe
        if not quiet:
            print(' normalmodes: running', exe, '...')
        if wiz is not None:
            wiz.message[1] = 'running diagrtb'
            cmd.refresh_wizard()
        process = subprocess.Popen([exe], cwd=tempdir,
                stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

        for line in process.stdout:
            if not quiet:
                sys.stdout.write(line)

        eigenfacs, frequencies = parse_eigenfacs(
                os.path.join(tempdir, 'diagrtb.eigenfacs'), last)

        if wiz is not None:
            wiz.message[1] = 'generating objects'
            cmd.refresh_wizard()

        states = int(states)
        factor = float(factor)
        if factor < 0:
            factor = natoms**0.5
        for mode in range(first, last+1):
            name = prefix + '%d' % mode
            cmd.delete(name)

            if not quiet:
                print(' normalmodes: object "%s" for mode %d with freq. %.6f' % \
                        (name, mode, frequencies[mode-1]))

            for state in range(1, states+1):
                cmd.create(name, sele_name, 1, state, zoom=0)
                cmd.alter_state(state, name,
                        '(x,y,z) = cpv.add([x,y,z], cpv.scale(next(myit), myfac))',
                        space={'cpv': cpv, 'myit': iter(eigenfacs[mode-1]),
                            'next': next,
                            'myfac': factor * (state - (states+1)/2.0)})

        # if CA only selection, show ribbon trace
        if natoms == cmd.count_atoms('(%s) and name CA' % sele_name):
            cmd.set('ribbon_trace_atoms', 1, prefix + '*')
            cmd.show_as('ribbon', prefix + '*')

        # store results
        if not hasattr(stored, 'nma_results'):
            stored.nma_results = []
        stored.nma_results.append({
            'facs': eigenfacs,
            'freq': frequencies,
            'sele': sele_name,
        })

    except OSError:
        print('Cannot execute "%s", please provide full path to executable' % (exe))
    except CmdException as e:
        print(' normalmodes: failed!', e)
    finally:
        if clean:
            shutil.rmtree(tempdir)
        elif not quiet:
            print(' normalmodes: Working directory "%s" not removed!' % (tempdir))
        # cmd.delete(sele_name)
        if wiz is not None:
            cmd.set_wizard_stack([w for w in cmd.get_wizard_stack() if w != wiz])

def parse_eigenfacs(filename='diagrtb.eigenfacs', readmax=20):
    line_it = iter(open(filename))
    eigenfacs = []
    values = []
    for line in line_it:
        a = line.split()
        if a[0] == 'VECTOR':
            assert a[2] == 'VALUE'
            number = int(a[1])
            assert number == len(eigenfacs) + 1
            if number > readmax:
                break
            value = float(a[3])
            next(line_it)
            vector = []
            eigenfacs.append(vector)
            values.append(value)
        else:
            a = list(map(float, a))
            vector.append(a)
    return eigenfacs, values

def normalmodes_mmtk(selection, cutoff=12.0, ff='Deformation', first=7, last=10,
        prefix='mmtk', states=7, factor=-1, quiet=1):
    '''
DESCRIPTION

    Fast normal modes for large proteins using an elastic network model (CA only)

    Based on:
    http://dirac.cnrs-orleans.fr/MMTK/using-mmtk/mmtk-example-scripts/normal-modes/
    '''
    try:
        import MMTK
    except ImportError:
        print('Failed to import MMTK, please add to PYTHONPATH')
        raise CmdException

    selection = selector.process(selection)
    cutoff = float(cutoff)
    first, last = int(first), int(last)
    states, factor, quiet = int(states), float(factor), int(quiet)

    from math import log
    from chempy import cpv

    from MMTK import InfiniteUniverse
    from MMTK.PDB import PDBConfiguration
    from MMTK.Proteins import Protein
    from MMTK.NormalModes import NormalModes

    from MMTK.ForceFields import DeformationForceField, CalphaForceField
    from MMTK.FourierBasis import FourierBasis, estimateCutoff
    from MMTK.NormalModes import NormalModes, SubspaceNormalModes

    model = 'calpha'
    ff = ff.lower()
    if 'deformationforcefield'.startswith(ff):
        forcefield = DeformationForceField(cutoff=cutoff/10.)
    elif 'calphaforcefield'.startswith(ff):
        forcefield = CalphaForceField(cutoff=cutoff/10.)
    elif 'amber94forcefield'.startswith(ff):
        from MMTK.ForceFields import Amber94ForceField
        forcefield = Amber94ForceField()
        model = 'all'
    else:
        raise NotImplementedError('unknown ff = ' + str(ff))
    if not quiet:
        print(' Forcefield:', forcefield.__class__.__name__)

    if model == 'calpha':
        selection = '(%s) and polymer and name CA' % (selection)

    f = StringIO(cmd.get_pdbstr(selection))
    conf = PDBConfiguration(f)
    items = conf.createPeptideChains(model)

    universe = InfiniteUniverse(forcefield)
    universe.protein = Protein(*items)

    nbasis = max(10, universe.numberOfAtoms()/5)
    cutoff, nbasis = estimateCutoff(universe, nbasis)
    if not quiet:
        print(" Calculating %d low-frequency modes." % nbasis)

    if cutoff is None:
        modes = NormalModes(universe)
    else:
        subspace = FourierBasis(universe, cutoff)
        modes = SubspaceNormalModes(universe, subspace)

    natoms = modes.array.shape[1]
    frequencies = modes.frequencies

    if factor < 0:
        factor = log(natoms)
        if not quiet:
            print(' set factor to %.2f' % (factor))

    if True: # cmd.count_atoms(selection) != natoms:
        import tempfile, os
        from MMTK import DCD
        filename = tempfile.mktemp(suffix='.pdb')
        sequence = DCD.writePDB(universe, None, filename)
        z = [a.index for a in sequence]
        selection = cmd.get_unused_name('_')
        cmd.load(filename, selection, zoom=0)
        os.remove(filename)

        if cmd.count_atoms(selection) != natoms:
            print('hmm... still wrong number of atoms')

    def eigenfacs_iter(mode):
        x = modes[mode-1].array
        return iter(x.take(z, 0))

    for mode in range(first, min(last, len(modes)) + 1):
        name = prefix + '%d' % mode
        cmd.delete(name)

        if not quiet:
            print(' normalmodes: object "%s" for mode %d with freq. %.6f' % \
                    (name, mode, frequencies[mode-1]))

        for state in range(1, states+1):
            cmd.create(name, selection, 1, state, zoom=0)
            cmd.alter_state(state, name,
                    '(x,y,z) = cpv.add([x,y,z], cpv.scale(next(myit), myfac))',
                    space={'cpv': cpv, 'myit': eigenfacs_iter(mode),
                        'next': next,
                        'myfac': 1e2 * factor * ((state-1.0)/(states-1.0) - 0.5)})

    cmd.delete(selection)
    if model == 'calpha':
        cmd.set('ribbon_trace_atoms', 1, prefix + '*')
        cmd.show_as('ribbon', prefix + '*')
    else:
        cmd.show_as('lines', prefix + '*')

def normalmodes_prody(selection, cutoff=15, first=7, last=10, guide=1,
        prefix='prody', states=7, factor=-1, quiet=1):
    '''
DESCRIPTION

    Anisotropic Network Model (ANM) analysis with ProDy.

    Based on:
    http://www.csb.pitt.edu/prody/examples/dynamics/enm/anm.html
    '''
    try:
        import prody
    except ImportError:
        print('Failed to import prody, please add to PYTHONPATH')
        raise CmdException

    first, last, guide = int(first), int(last), int(guide)
    states, factor, quiet = int(states), float(factor), int(quiet)
    assert first > 6

    if guide:
        selection = '(%s) and guide and alt A+' % (selection)
    tmpsele = cmd.get_unused_name('_')
    cmd.select(tmpsele, selection)

    f = StringIO(cmd.get_pdbstr(tmpsele))
    conf = prody.parsePDBStream(f)

    modes = prody.ANM()
    modes.buildHessian(conf, float(cutoff))
    modes.calcModes(last - first + 1)

    if factor < 0:
        from math import log
        natoms = modes.numAtoms()
        factor = log(natoms) * 10
        if not quiet:
            print(' set factor to %.2f' % (factor))

    for mode in range(first, last + 1):
        name = prefix + '%d' % mode
        cmd.delete(name)

        if not quiet:
            print(' normalmodes: object "%s" for mode %d' % (name, mode))

        for state in range(1, states+1):
            xyz_it = iter(modes[mode-7].getArrayNx3() * (factor *
                    ((state-1.0)/(states-1.0) - 0.5)))
            cmd.create(name, tmpsele, 1, state, zoom=0)
            cmd.alter_state(state, name, '(x,y,z) = next(xyz_it) + (x,y,z)',
                    space={'xyz_it': xyz_it, 'next': next})

    cmd.delete(tmpsele)

    if guide:
        cmd.set('ribbon_trace_atoms', 1, prefix + '*')
        cmd.show_as('ribbon', prefix + '*')
    else:
        cmd.show_as('lines', prefix + '*')

cmd.extend('normalmodes_pdbmat', normalmodes_pdbmat)
cmd.extend('normalmodes_mmtk', normalmodes_mmtk)
cmd.extend('normalmodes_prody', normalmodes_prody)

cmd.auto_arg[0].update([
    ('normalmodes_pdbmat', cmd.auto_arg[0]['zoom']),
    ('normalmodes_mmtk', cmd.auto_arg[0]['zoom']),
    ('normalmodes_prody', cmd.auto_arg[0]['zoom']),
])

# vi: ts=4:sw=4:smarttab:expandtab
