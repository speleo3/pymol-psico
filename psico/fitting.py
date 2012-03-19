'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def save_pdb_without_ter(filename, selection, **kwargs):
    '''
DESCRIPTION

    Save PDB file without TER records. External applications like TMalign and
    DynDom stop reading PDB files at TER records, which might be undesired in
    case of missing loops.
    '''
    v = cmd.get_setting_boolean('pdb_use_ter_records')
    if v: cmd.unset('pdb_use_ter_records')
    cmd.save(filename, selection, **kwargs)
    if v: cmd.set('pdb_use_ter_records')

def alignwithanymethod(mobile, target, methods=None, async=1, quiet=1):
    '''
DESCRIPTION

    Align copies of mobile to target with several alignment methods

ARGUMENTS

    mobile = string: atom selection

    target = string: atom selection

    methods = string: space separated list of PyMOL commands which take
    arguments "mobile" and "target" (in any order) {default: align super
    cealign bfit tmalign}
    '''
    import threading
    import time
    if methods is None:
        methods = align_methods
    else:
        methods = methods.split()
    async, quiet = int(async), int(quiet)
    mobile_obj = cmd.get_object_list('first (' + mobile + ')')[0]
    def myalign(method):
        newmobile = cmd.get_unused_name(mobile_obj + '_' + method)
        cmd.create(newmobile, mobile_obj)
        start = time.time()
        cmd.do('%s mobile=%s in %s, target=%s' % (method, newmobile, mobile, target))
        if not quiet:
            print 'Finished: %s (%.2f sec)' % (method, time.time() - start)
    for method in methods:
        if method not in cmd.keyword:
            if not quiet:
                print 'No such method:', method
            continue
        if async:
            t = threading.Thread(target=myalign, args=(method,))
            t.setDaemon(1)
            t.start()
        else:
            myalign(method)

def tmalign(mobile, target, mobile_state=1, target_state=1, args='',
        exe='TMalign', ter=0, transform=1, object=None, quiet=0):
    '''
DESCRIPTION

    TMalign wrapper

    Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
    http://zhanglab.ccmb.med.umich.edu/TM-align/

ARGUMENTS

    mobile, target = string: atom selections

    mobile_state, target_state = int: object states {default: 1}

    args = string: Extra arguments like -d0 5 -L 100

    exe = string: Path to TMalign executable {default: TMalign}

    ter = 0/1: If ter=0, then ignore chain breaks because TMalign will stop
    at first TER record {default: 0}

SEE ALSO

    tmscore, mmalign
    '''
    import subprocess, tempfile, os, re

    ter, quiet = int(ter), int(quiet)

    mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    target_filename = tempfile.mktemp('.pdb', 'target')
    mobile_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (mobile)
    target_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (target)

    if ter:
        save = cmd.save
    else:
        save = save_pdb_without_ter
    save(mobile_filename, mobile_ca_sele, state=mobile_state)
    save(target_filename, target_ca_sele, state=target_state)

    exe = cmd.exp_path(exe)
    args = [exe, mobile_filename, target_filename] + args.split()

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE)
        lines = process.stdout.readlines()
    except OSError:
        print 'Cannot execute "%s", please provide full path to TMscore or TMalign executable' % (exe)
        raise CmdException
    finally:
        os.remove(mobile_filename)
        os.remove(target_filename)

    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    matrix = []
    line_it = iter(lines)
    headercheck = False
    for line in line_it:
        if 4 >= rowcount > 0:
            if rowcount >= 2:
                a = map(float, line.split())
                matrix.extend(a[2:5])
                matrix.append(a[1])
            rowcount += 1
        elif not headercheck and line.startswith(' * '):
            a = line.split(None, 2)
            if len(a) == 3:
                headercheck = a[1]
        elif line.lower().startswith(' -------- rotation matrix'):
            rowcount = 1
        elif line.startswith('(":" denotes'):
            alignment = [line_it.next().rstrip() for i in range(3)]
        else:
            match = re_score.search(line)
            if match is not None:
                r = float(match.group(1))
        if not quiet:
            print line.rstrip()
        
    if not quiet:
        for i in range(0, len(alignment[0])-1, 78):
            for line in alignment:
                print line[i:i+78]
            print ''

    assert len(matrix) == 3*4
    matrix.extend([0,0,0,1])

    if int(transform):
        for model in cmd.get_object_list('(' + mobile + ')'):
            cmd.transform_object(model, matrix, state=0, homogenous=1)
    
    # alignment object
    if object is not None:
        mobile_idx, target_idx = [], []
        space = {'mobile_idx': mobile_idx, 'target_idx': target_idx}
        cmd.iterate(mobile_ca_sele, 'mobile_idx.append("%s`%d" % (model, index))', space=space)
        cmd.iterate(target_ca_sele, 'target_idx.append("%s`%d" % (model, index))', space=space)
        for i, aa in enumerate(alignment[0]):
            if aa == '-':
                mobile_idx.insert(i, None)
        for i, aa in enumerate(alignment[2]):
            if aa == '-':
                target_idx.insert(i, None)
        if (len(mobile_idx) == len(target_idx) == len(alignment[2])):
            cmd.rms_cur(
                    ' '.join(idx for (idx, m) in zip(mobile_idx, alignment[1]) if m in ':.'),
                    ' '.join(idx for (idx, m) in zip(target_idx, alignment[1]) if m in ':.'),
                    cycles=0, matchmaker=4, object=object)
        else:
            print 'Could not load alignment object'

    if not quiet:
        if headercheck:
            print 'Finished Program:', headercheck
        if r is not None:
            print 'Found in output TM-score = %.4f' % (r)

    return r

def tmscore(mobile, target, mobile_state=1, target_state=1, args='',
        exe='TMscore', ter=0, transform=1, object=None, quiet=0):
    '''
DESCRIPTION

    TMscore wrapper

    Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710
    http://zhanglab.ccmb.med.umich.edu/TM-score/

ARGUMENTS

    mobile, target = string: atom selections

    mobile_state, target_state = int: object states {default: 1}

    args = string: Extra arguments like -d 5

    exe = string: Path to TMscore executable {default: TMscore}

    ter = 0/1: If ter=0, then ignore chain breaks because TMscore will stop
    at first TER record {default: 0}

SEE ALSO

    tmalign, mmalign
    '''
    return tmalign(mobile, target, mobile_state, target_state, args,
            exe, ter, transform, object, quiet)

def mmalign(mobile, target, mobile_state=1, target_state=1, args='',
        exe='MMalign', ter=0, transform=1, quiet=0):
    '''
DESCRIPTION

    MMalign wrapper

    Reference: S. Mukherjee and Y. Zhang, Nucleic Acids Research 2009; 37: e83
    http://zhanglab.ccmb.med.umich.edu/MM-align/

SEE ALSO

    tmalign, tmscore
    '''
    return tmalign(mobile, target, mobile_state, target_state, args,
            exe, ter, transform, None, quiet)

def dyndom_parse_info(filename, selection='(all)', quiet=0):
    import re
    fixed = False
    fixed_name = None
    dom_nr = 0
    color = 'none'
    bending = list()
    for line in open(filename):
        if line.startswith('FIXED  DOMAIN'):
            fixed = True
            continue
        if line.startswith('MOVING DOMAIN'):
            fixed = False
            continue
        m = re.match(r'DOMAIN NUMBER: *(\d+) \(coloured (\w+)', line)
        if m:
            dom_nr = m.group(1)
            color = m.group(2)
            continue
        m = re.match(r'RESIDUE NUMBERS :(.*)', line)
        if m:
            resi = m.group(1)
            resi = resi.replace(',', '+')
            resi = resi.replace(' ', '')
            if not quiet:
                print 'Domain ' + dom_nr + ' (' + color + '): resi ' + resi
            name = 'domain_' + dom_nr
            cmd.select(name, '(%s) and (resi %s)' % (selection, resi))
            cmd.color(color, name)
            if fixed:
                fixed_name = name
            continue
        m = re.match(r'BENDING RESIDUES:(.*)', line)
        if m:
            resi = m.group(1)
            resi = resi.replace(',', '+')
            resi = resi.replace(' ', '')
            bending.append(resi)
    if len(bending) > 0:
        name = 'bending'
        cmd.select(name, '(%s) and (resi %s)' % (selection, '+'.join(bending)))
        cmd.color('green', name)
    return fixed_name

def dyndom(mobile, target, window=5, domain=20, ratio=1.0, exe='DynDom', transform=1, quiet=1):
    '''
DESCRIPTION

    DynDom wrapper

    DynDom is a program to determine domains, hinge axes and hinge bending
    residues in proteins where two conformations are available.

    http://fizz.cmp.uea.ac.uk/dyndom/

USAGE

    dyndom mobile, target [, window [, domain [, ratio ]]]

SEE ALSO

    promix
    '''
    import tempfile, subprocess, os, shutil, sys

    window, domain, ratio = int(window), int(domain), float(ratio)
    transform, quiet = int(transform), int(quiet)

    chains = cmd.get_chains(mobile)
    if len(chains) != 1:
        print 'mobile selection must be single chain'
        raise CmdException
    chain1id = chains[0]
    chains = cmd.get_chains(target)
    if len(chains) != 1:
        print 'target selection must be single chain'
        raise CmdException
    chain2id = chains[0]

    exe = cmd.exp_path(exe)
    tempdir = tempfile.mkdtemp()

    try:
        filename1 = os.path.join(tempdir, 'mobile.pdb')
        filename2 = os.path.join(tempdir, 'target.pdb')
        commandfile = os.path.join(tempdir, 'command.txt')
        infofile = os.path.join(tempdir, 'out_info')

        save_pdb_without_ter(filename1, '(%s) and polymer' % mobile)
        save_pdb_without_ter(filename2, '(%s) and polymer' % target)

        f = open(commandfile, 'w')
        f.write('title=out\nfilename1=%s\nchain1id=%s\nfilename2=%s\nchain2id=%s\n' \
                'window=%d\ndomain=%d\nratio=%.4f\n' % (filename1, chain1id,
                    filename2, chain2id, window, domain, ratio))
        f.close()

        process = subprocess.Popen([exe, commandfile], cwd=tempdir,
                stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

        for line in process.stdout:
            if not quiet:
                sys.stdout.write(line)

        cmd.color('gray', mobile)
        fixed_name = dyndom_parse_info(infofile, mobile, quiet)
    except OSError:
        print 'Cannot execute "%s", please provide full path to DynDom executable' % (exe)
        raise CmdException
    finally:
        shutil.rmtree(tempdir)

    if transform and fixed_name is not None:
        cmd.align(fixed_name, target)

def gdt_ts(mobile, target, cutoffs='1 2 4 8', quiet=1):
    '''
DESCRIPTION

    Global Distance Test Total Score (GDT_TS)
    '''
    cutoffs = map(float, cutoffs.split())
    quiet = int(quiet)
    mobile = '(' + mobile + ') and guide'
    target = '(' + target + ') and guide'
    ts = 0
    N = min(cmd.count_atoms(mobile), cmd.count_atoms(target))
    for cutoff in cutoffs:
        x = cmd.align(mobile, target, cutoff=cutoff, transform=0)
        p = float(x[1]) / N
        if not quiet:
            print ' GDT_TS: GDT_P%.1f = %.2f' % (cutoff, p)
        ts += p
    ts /= len(cutoffs)
    if not quiet:
        print ' GDT_TS: Total Score = %.2f' % (ts)
    return ts

def get_rmsd_func():
    '''
DESCRIPTION

    API only. Returns a function that uses either numpy (fast) or chempy.cpv
    (slow) to calculate the rmsd fit of two nx3 arrays.
    '''
    try:
        # this is much faster than cpv.fit
        from numpy import dot, sqrt, array
        from numpy.linalg import svd
        def rmsd(X, Y):
            X = X - X.mean(0)
            Y = Y - Y.mean(0)
            R_x = (X**2).sum()
            R_y = (Y**2).sum()
            L = svd(dot(Y.T, X))[1]
            return sqrt((R_x + R_y - 2 * L.sum()) / len(X))
        rmsd.array = array
    except ImportError:
        from chempy import cpv
        def rmsd(X, Y):
            return cpv.fit(X, Y)[-1]
        rmsd.array = lambda x: x
    return rmsd

def matchmaker(mobile, target, match):
    '''
DESCRIPTION

    API only. See "local_rms" for description of "match" argument.
    '''
    temporary_names = []
    if match == 'none':
        pass
    elif match in ['in', 'like']:
        mobile = '(%s) %s (%s)' % (mobile, match, target)
        target = '(%s) %s (%s)' % (target, match, mobile)
    elif match in ['align', 'super']:
        aln_obj = cmd.get_unused_name('_')
        temporary_names.append(aln_obj)
        align = cmd.keyword[match][0]
        align(mobile, target, cycles=0, transform=0, object=aln_obj)
        from .selecting import wait_for
        wait_for(aln_obj)
        mobile = '(%s) and %s' % (mobile, aln_obj)
        target = '(%s) and %s' % (target, aln_obj)
    elif match in cmd.get_names('all') and cmd.get_type(match) == 'object:':
        names = cmd.get_object_list('(' + match + ')')
        if len(names) == 2:
            # easy case, alignment object covers only the two given objects
            mobile = '(%s) and %s' % (mobile, match)
            target = '(%s) and %s' % (target, match)
        else:
            # if alignment object covers more than the two objects, then we
            # need to pick those columns that have no gap in any of the two
            # given objects
            mobile_object = cmd.get_object_list('(' + mobile + ')')[0]
            target_object = cmd.get_object_list('(' + target + ')')[0]
            mobile_index_list = []
            target_index_list = []
            for column in cmd.get_raw_alignment(match):
                column = dict(column)
                if mobile_object in column and target_object in column:
                    mobile_index_list.append(column[mobile_object])
                    target_index_list.append(column[target_object])
            mobile_name = cmd.get_unused_name('_mobile')
            target_name = cmd.get_unused_name('_target')
            cmd.select_list(mobile_name, mobile_object, mobile_index_list, mode='index')
            cmd.select_list(target_name, target_object, target_index_list, mode='index')
            temporary_names.append(mobile_name)
            temporary_names.append(target_name)
            mobile = '(%s) and %s' % (mobile, mobile_name)
            target = '(%s) and %s' % (target, target_name)
    else:
        print ' Error: unkown match method', match
        raise CmdException
    return mobile, target, temporary_names

def local_rms(mobile, target, window=20, mobile_state=1, target_state=1,
        match='align', load_b=1, visualize=1, quiet=1):
    '''
DESCRIPTION

    "local_rms" computes the C-alpha RMS fit within a sliding window along the
    backbone. The obtained RMS is assigned as a pseudo b-factor to the residue
    in the middle of the window. This is useful to visualize hinge-regions.

    The result is very sensitive to window size.

USAGE

    local_rms mobile, target [, window ]

ARGUMENTS

    mobile = string: object to assign b-factors and to visualize as putty cartoon

    target = string: object to superimpose mobile to

    window = integer: width of sliding window {default: 20}

    match = string: in, like, align, none or the name of an alignment object
    {default: align}
      * in: match all atom identifiers (segi,chain,resn,resi,name)
      * like: match residue number (resi)
      * align: do a sequence alignment
      * none: assume same number of atoms in both selections
      * name of alignment object: take sequence alignment from object

EXAMPLE

    fetch 2x19 2xwu, async=0
    remove not chain B or not polymer
    local_rms 2x19, 2xwu, 40
    '''
    rmsd = get_rmsd_func()
    array = rmsd.array

    window = int(window)
    mobile_state, target_state = int(mobile_state), int(target_state)
    load_b, visualize, quiet = int(load_b), int(visualize), int(quiet)

    w2 = window/2
    w4 = window/4

    mguide, tguide, tmp_names = matchmaker('(%s) and guide' % (mobile),
            '(%s) and guide' % (target), match)

    model_mobile = cmd.get_model(mguide)
    model_target = cmd.get_model(tguide)

    for name in tmp_names:
        cmd.delete(name)

    if len(model_mobile.atom) != len(model_mobile.atom):
        print 'Error: number of atoms differ, please check match method'
        raise CmdException

    seq_start = model_mobile.atom[0].resi_number
    seq_end = model_mobile.atom[-1].resi_number

    resv2i = dict((a.resi_number,i) for (i,a) in enumerate(model_mobile.atom))
    resv2b = dict()

    X_mobile = array(model_mobile.get_coord_list())
    X_target = array(model_target.get_coord_list())

    for resv in range(seq_start, seq_end + 1):
        for resv_from in range(resv-w2, resv+1):
            i_from = resv2i.get(resv_from)
            if i_from is not None:
                break
        for resv_to in range(resv+w2, resv-1, -1):
            i_to = resv2i.get(resv_to)
            if i_to is not None:
                break
        if i_from is None or i_to is None:
            continue
        if i_to - i_from < w4:
            continue

        x = X_mobile[i_from:i_to+1]
        y = X_target[i_from:i_to+1]
        resv2b[resv] = rmsd(x, y)

        if not quiet:
            print ' resi %4d: RMS = %6.3f (%4d atoms)' % (resv, resv2b[resv], i_to - i_from + 1)

    if load_b:
        cmd.alter(mobile, 'b=resv2b.get(resv, -1.0)', space={'resv2b': resv2b})

    if load_b and visualize:
        cmd.color('yellow', '(%s) and b < -0.5' % (mobile))
        cmd.spectrum('b', 'blue_white_red', '(%s) and b > -0.5' % (mobile))
        cmd.show_as('cartoon', mobile)
        cmd.hide('cartoon', '(%s) and b < -0.5' % (mobile))
        cmd.cartoon('putty', mobile)

    return resv2b

def extra_fit(selection='(all)', reference=None, method='align', zoom=1,
        quiet=0, _self=cmd, **kwargs):
    '''
DESCRIPTION

    Like "intra_fit", but for multiple objects instead of
    multiple states.

ARGUMENTS

    selection = string: atom selection of multiple objects {default: all}

    reference = string: reference object name {default: first object in selection}

    method = string: alignment method (command that takes "mobile" and "target"
    arguments, like "align", "super", "cealign" {default: align}

    ... extra arguments are passed to "method"

SEE ALSO

    alignto, cmd.util.mass_align, align_all.py from Robert Campbell
    '''
    zoom, quiet = int(zoom), int(quiet)
    sele_name = cmd.get_unused_name('_')
    cmd.select(sele_name, selection) # for speed
    models = cmd.get_object_list(sele_name)
    if reference is None:
        reference = models[0]
        models = models[1:]
    elif reference in models:
        models.remove(reference)
    else:
        cmd.select(sele_name, reference, merge=1)
    if cmd.is_string(method):
        if method in cmd.keyword:
            method = cmd.keyword[method][0]
        else:
            print 'Unknown method:', method
            raise CmdException
    for model in models:
        x = method(mobile='%s and model %s' % (sele_name, model),
                target='%s and model %s' % (sele_name, reference), **kwargs)
        if not quiet:
            if cmd.is_sequence(x):
                print '%-20s RMS = %8.3f (%d atoms)' % (model, x[0], x[1])
            elif isinstance(x, float):
                print '%-20s RMS = %8.3f' % (model, x)
    if zoom:
        cmd.zoom(sele_name)
    cmd.delete(sele_name)

def _run_theseus(args, tempdir, preserve, quiet):
    '''
DESCRIPTION

    Helper function for theseus and intra_theseus
    '''
    import subprocess, os

    translations = []
    rotations = []

    try:
        if quiet:
            subprocess.call(args, cwd=tempdir)
        else:
            import re
            unesc = re.compile('\x1b' + r'\[[\d;]+m').sub

            process = subprocess.Popen(args, cwd=tempdir, stdout=subprocess.PIPE)
            for line in process.stdout:
                print unesc('', line.rstrip())

        handle = open(os.path.join(tempdir, 'theseus_transf2.txt'))
        for line in handle:
            if line[10:13] == ' t:':
                translations.append(map(float, line[13:].split()))
            elif line[10:13] == ' R:':
                rotations.append(map(float, line[13:].split()))
        handle.close()

    except OSError:
        print ' Error: Cannot execute "%s"' % (args[0])
        raise CmdException
    finally:
        if not preserve:
            import shutil
            shutil.rmtree(tempdir)
        elif not quiet:
            print ' Not deleting temporary directory:', tempdir

    return translations, rotations

def theseus(mobile, target, match='align', cov=0, cycles=200,
        mobile_state=1, target_state=1, exe='theseus', preserve=0, quiet=1):
    '''
DESCRIPTION

    Structural superposition of two molecules with maximum likelihood.

    THESEUS: Maximum likelihood multiple superpositioning
    http://www.theseus3d.org

ARGUMENTS

    mobile = string: atom selection for mobile atoms

    target = string: atom selection for target atoms

    match = string: in, like, align, none or the name of an alignment object
    (see "local_rms" help for details) {default: align}

    cov = 0/1: 0 is variance weighting, 1 is covariance weighting (slower)
    {default: 0}

SEE ALSO

    align, super, cealign
    '''
    import tempfile, os

    cov, cycles = int(cov), int(cycles)
    mobile_state, target_state = int(mobile_state), int(target_state)
    preserve, quiet = int(preserve), int(quiet)

    tempdir = tempfile.mkdtemp()
    mobile_filename = os.path.join(tempdir, 'mobile.pdb')
    target_filename = os.path.join(tempdir, 'target.pdb')

    mmobile, mtarget, tmp_names = matchmaker(mobile, target, match)
    cmd.save(mobile_filename, mmobile, mobile_state)
    cmd.save(target_filename, mtarget, target_state)
    for name in tmp_names:
        cmd.delete(name)

    exe = cmd.exp_path(exe)
    args = [exe, '-a0', '-c' if cov else '-v', '-i%d' % cycles,
            mobile_filename, target_filename]

    translations, rotations = _run_theseus(args, tempdir, preserve, quiet)
    matrices = [R[0:3] + [i*t[0]] + R[3:6] + [i*t[1]] + R[6:9] + [i*t[2], 0,0,0, 1]
            for (R, t, i) in zip(rotations, translations, [-1,1])]

    obj_list = cmd.get_object_list('(' + mobile + ')')
    for obj in obj_list:
        cmd.transform_object(obj, matrices[0], 0, transpose=1)
        cmd.transform_object(obj, matrices[1], 0)

    if not quiet:
        print ' theseus: done'

def intra_theseus(selection, state=1, cov=0, cycles=200,
        exe='theseus', preserve=0, quiet=1):
    '''
DESCRIPTION

    Fits all states of an object to an atom selection with maximum likelihood.

    THESEUS: Maximum likelihood multiple superpositioning
    http://www.theseus3d.org

ARGUMENTS

    selection = string: atoms to fit

    state = integer: keep transformation of this state unchanged {default: 1}

    cov = 0/1: 0 is variance weighting, 1 is covariance weighting (slower)
    {default: 0}

SEE ALSO

    intra_fit, intra_rms_cur
    '''
    import tempfile, os

    state, cov, cycles = int(state), int(cov), int(cycles)
    preserve, quiet = int(preserve), int(quiet)

    tempdir = tempfile.mkdtemp()
    filename = os.path.join(tempdir, 'mobile.pdb')

    cmd.save(filename, selection, 0)

    exe = cmd.exp_path(exe)
    args = [exe, '-a0', '-c' if cov else '-v', '-i%d' % cycles, filename]

    translations = []
    rotations = []

    translations, rotations = _run_theseus(args, tempdir, preserve, quiet)
    matrices = [R[0:3] + [-t[0]] + R[3:6] + [-t[1]] + R[6:9] + [-t[2], 0,0,0, 1]
            for (R, t) in zip(rotations, translations)]

    # intra fit states
    obj_list = cmd.get_object_list('(' + selection + ')')
    for i, m in enumerate(matrices):
        for obj in obj_list:
            cmd.transform_object(obj, m, i+1, transpose=1)

    # fit back to given state
    if 0 < state <= len(matrices):
        m = list(matrices[state-1])
        for i in [3,7,11]:
            m[i] *= -1
        for obj in obj_list:
            cmd.transform_object(obj, m, 0)

    if not quiet:
        print ' intra_theseus: %d states aligned' % (len(matrices))

# all those have kwargs: mobile, target, mobile_state, target_state
align_methods = ['align', 'super', 'cealign', 'tmalign', 'theseus']
align_methods_sc = cmd.Shortcut(align_methods)

# pymol commands
cmd.extend('alignwithanymethod', alignwithanymethod)
cmd.extend('tmalign', tmalign)
cmd.extend('tmscore', tmscore)
cmd.extend('mmalign', mmalign)
cmd.extend('dyndom', dyndom)
cmd.extend('gdt_ts', gdt_ts)
cmd.extend('local_rms', local_rms)
cmd.extend('extra_fit', extra_fit)
cmd.extend('intra_theseus', intra_theseus)
cmd.extend('theseus', theseus)

# autocompletion
cmd.auto_arg[0].update({
    'alignwithanymethod': cmd.auto_arg[0]['align'],
    'tmalign': cmd.auto_arg[0]['align'],
    'tmscore': cmd.auto_arg[0]['align'],
    'mmalign': cmd.auto_arg[0]['align'],
    'dyndom': cmd.auto_arg[0]['align'],
    'gdt_ts': cmd.auto_arg[0]['align'],
    'local_rms': cmd.auto_arg[0]['align'],
    'extra_fit': cmd.auto_arg[0]['align'],
    'theseus': cmd.auto_arg[0]['align'],
    'intra_theseus': cmd.auto_arg[0]['align'],
})
cmd.auto_arg[1].update({
    'alignwithanymethod': cmd.auto_arg[1]['align'],
    'tmalign': cmd.auto_arg[1]['align'],
    'tmscore': cmd.auto_arg[1]['align'],
    'mmalign': cmd.auto_arg[1]['align'],
    'dyndom': cmd.auto_arg[1]['align'],
    'gdt_ts': cmd.auto_arg[1]['align'],
    'local_rms': cmd.auto_arg[1]['align'],
    'extra_fit': cmd.auto_arg[0]['disable'],
    'theseus': cmd.auto_arg[1]['align'],
})
cmd.auto_arg[2].update([
    ('extra_fit', [ align_methods_sc, 'alignment method', '' ]),
])

# vi: ts=4:sw=4:smarttab:expandtab
