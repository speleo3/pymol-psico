'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def alignwithanymethod(mobile, target, methods=None, async=1, quiet=1):
    '''
DESCRIPTION

    Align copies of mobile to target with several alignment methods

ARGUMENTS

    mobile = string: atom selection

    target = string: atom selection

    methods = string: space separated list of PyMOL commands which take
    arguments "mobile" and "target" (in any order) {default: align super
    cealign tmalign theseus}
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

    TMalign wrapper. You may also use this as a TMscore or MMalign wrapper
    if you privide the corresponding executable with the "exe" argument.

    Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
    http://zhanglab.ccmb.med.umich.edu/TM-align/

ARGUMENTS

    mobile, target = string: atom selections

    mobile_state, target_state = int: object states {default: 1}

    args = string: Extra arguments like -d0 5 -L 100

    exe = string: Path to TMalign (or TMscore, MMalign) executable
    {default: TMalign}

    ter = 0/1: If ter=0, then ignore chain breaks because TMalign will stop
    at first TER record {default: 0}
    '''
    import subprocess, tempfile, os, re
    from .exporting import save_pdb_without_ter

    ter, quiet = int(ter), int(quiet)

    mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    target_filename = tempfile.mktemp('.pdb', 'target')
    matrix_filename = tempfile.mktemp('.txt', 'matrix')
    mobile_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (mobile)
    target_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (target)

    if ter:
        save = cmd.save
    else:
        save = save_pdb_without_ter
    save(mobile_filename, mobile_ca_sele, state=mobile_state)
    save(target_filename, target_ca_sele, state=target_state)

    exe = cmd.exp_path(exe)
    args = [exe, mobile_filename, target_filename, '-m', matrix_filename] + args.split()

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE)
        lines = process.stdout.readlines()
    except OSError:
        print 'Cannot execute "%s", please provide full path to TMscore or TMalign executable' % (exe)
        raise CmdException
    finally:
        os.remove(mobile_filename)
        os.remove(target_filename)

    # TMalign >= 2012/04/17
    if os.path.exists(matrix_filename):
        lines += open(matrix_filename).readlines()
        os.remove(matrix_filename)

    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    matrix = []
    line_it = iter(lines)
    headercheck = False
    alignment = []
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

def dyndom(mobile, target, window=5, domain=20, ratio=1.0, exe='', transform=1, quiet=1):
    '''
DESCRIPTION

    DynDom wrapper

    DynDom is a program to determine domains, hinge axes and hinge bending
    residues in proteins where two conformations are available.

    http://fizz.cmp.uea.ac.uk/dyndom/

USAGE

    dyndom mobile, target [, window [, domain [, ratio ]]]
    '''
    import tempfile, subprocess, os, shutil, sys
    from .exporting import save_pdb_without_ter

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

    if not exe:
        from . import which
        exe = which('DynDom', 'dyndom')
        if not exe:
            print ' Error: Cannot find DynDom executable'
            raise CmdException
    else:
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

    Legacy, deprecated, use MatchMaker instead
    '''
    mm = MatchMaker(mobile, target, match)
    mm.autodelete = False
    return mm.mobile, mm.target, mm.temporary

class MatchMaker(object):
    '''
DESCRIPTION

    API only. Matches two atom selections and provides two matched
    subselections with equal atom count. May involve temporary objects
    or named selections which will be automatically deleted.

ARGUMENTS

    mobile = string: first atom selection

    target = string: second atom selection

    match = string: method how to match atoms
        * none: (dummy)
        * in: match atoms by "in" operator
        * like: match atoms by "like" operator
        * align: match atoms by cmd.align (without refinement)
        * super: match atoms by cmd.super (without refinement)
        * <name of alignment object>: use given alignment

RESULT

    Properties "mobile" and "target" hold the matched subselections as
    selection strings.
    '''
    def __init__(self, mobile, target, match):
        self.autodelete = True
        self.temporary = []

        if match == 'none':
            self.mobile = mobile
            self.target = target
        elif match in ['in', 'like']:
            self.mobile = '(%s) %s (%s)' % (mobile, match, target)
            self.target = '(%s) %s (%s)' % (target, match, mobile)
        elif match in ['align', 'super']:
            self.align(mobile, target, match)
        elif match in cmd.get_names('all') and cmd.get_type(match) in ('object:', 'object:alignment'):
            self.from_alignment(mobile, target, match)
        else:
            print ' Error: unkown match method', match
            raise CmdException

    def check(self):
        return cmd.count_atoms(self.mobile) == cmd.count_atoms(self.target)

    def align(self, mobile, target, match):
        '''
        Align mobile to target using the alignment method given by "match"
        '''
        aln_obj = cmd.get_unused_name('_')
        self.temporary.append(aln_obj)

        align = cmd.keyword[match][0]
        align(mobile, target, cycles=0, transform=0, object=aln_obj)
        cmd.disable(aln_obj)

        self.from_alignment(mobile, target, aln_obj)

    def from_alignment(self, mobile, target, aln_obj):
        '''
        Use alignment given by "aln_obj" (name of alignment object)
        '''
        from .selecting import wait_for
        wait_for(aln_obj)

        self.mobile = '(%s) and %s' % (mobile, aln_obj)
        self.target = '(%s) and %s' % (target, aln_obj)
        if self.check():
            return

        # difficult: if selections spans only part of the alignment or
        # if alignment object covers more than the two objects, then we
        # need to pick those columns that have no gap in any of the two
        # given selections

        mobileidx = set(cmd.index(mobile))
        targetidx = set(cmd.index(target))
        mobileidxsel = []
        targetidxsel = []

        for column in cmd.get_raw_alignment(aln_obj):
            mobiles = mobileidx.intersection(column)
            if len(mobiles) == 1:
                targets = targetidx.intersection(column)
                if len(targets) == 1:
                    mobileidxsel.extend(mobiles)
                    targetidxsel.extend(targets)

        self.mobile = cmd.get_unused_name('_mobile')
        self.target = cmd.get_unused_name('_target')
        self.temporary.append(self.mobile)
        self.temporary.append(self.target)

        mobile_objects = set(idx[0] for idx in mobileidxsel)
        target_objects = set(idx[0] for idx in targetidxsel)

        if len(mobile_objects) == len(target_objects) == 1:
            mobile_index_list = [idx[1] for idx in mobileidxsel]
            target_index_list = [idx[1] for idx in targetidxsel]
            cmd.select_list(self.mobile, mobile_objects.pop(), mobile_index_list, mode='index')
            cmd.select_list(self.target, target_objects.pop(), target_index_list, mode='index')
        else:
            cmd.select(self.mobile, ' '.join('%s`%d' % idx for idx in mobileidxsel))
            cmd.select(self.target, ' '.join('%s`%d' % idx for idx in targetidxsel))

    def __del__(self):
        if not self.autodelete:
            return
        for name in self.temporary:
            cmd.delete(name)

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

    mm = MatchMaker('(%s) and guide' % (mobile),
            '(%s) and guide' % (target), match)

    model_mobile = cmd.get_model(mm.mobile)
    model_target = cmd.get_model(mm.target)

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
    t_type = float

    try:
        if quiet:
            subprocess.call(args, cwd=tempdir)
        else:
            import re
            unesc = re.compile('\x1b' + r'\[[\d;]+m').sub

            process = subprocess.Popen(args, cwd=tempdir, stdout=subprocess.PIPE)
            for line in process.stdout:
                print unesc('', line.rstrip())

        filename = os.path.join(tempdir, 'theseus_transf2.txt')

        if not os.path.exists(filename):
            # THESEUS 3.x
            filename = os.path.join(tempdir, 'theseus_transf.txt')
            if not os.path.exists(filename):
                raise CmdException('no theseus_transf2.txt or '
                        'theseus_transf.txt output file')
            t_type = lambda t: float(t) * -1.

        handle = open(filename)
        for line in handle:
            if line[10:13] == ' t:':
                translations.append(map(t_type, line[13:].split()))
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

    mm = MatchMaker(mobile, target, match)
    cmd.save(mobile_filename, mm.mobile, mobile_state)
    cmd.save(target_filename, mm.target, target_state)

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

def prosmart(mobile, target, mobile_state=1, target_state=1,
        exe='prosmart', transform=1, object=None, quiet=0):
    '''
DESCRIPTION

    ProSMART wrapper.

    http://www2.mrc-lmb.cam.ac.uk/groups/murshudov/
    '''
    import subprocess, tempfile, os, shutil, glob

    quiet = int(quiet)

    tempdir = tempfile.mkdtemp()
    mobile_filename = os.path.join(tempdir, 'mobile.pdb')
    target_filename = os.path.join(tempdir, 'target.pdb')

    cmd.save(mobile_filename, mobile, state=mobile_state)
    cmd.save(target_filename, target, state=target_state)

    exe = cmd.exp_path(exe)
    args = [exe, '-p1', mobile_filename, '-p2', target_filename, '-a']

    xglob = lambda x: glob.glob(os.path.join(tempdir, 'ProSMART_Output/Output_Files', x))

    try:
        subprocess.check_call(args, cwd=tempdir)

        transfiles = xglob('Superposition/Transformations/*/*.txt')
        with open(transfiles[0]) as f:
            f = iter(f)
            for line in f:
                if line.startswith('ROTATION'):
                    matrix = [map(float, next(f).split()) + [0] for _ in range(3)]
                elif line.startswith('TRANSLATION'):
                    matrix.append([-float(v) for v in next(f).split()] + [1])
                    break

        if int(transform):
            matrix = [v for m in matrix for v in m]
            assert len(matrix) == 4*4
            for model in cmd.get_object_list('(' + mobile + ')'):
                cmd.transform_object(model, matrix, state=0)

        if object:
            from .importing import load_aln
            alnfiles = xglob('Residue_Alignment_Scores/*/*.txt')
            alnfiles = filter(lambda x: not x.endswith('_clusters.txt'), alnfiles)
            load_aln(alnfiles[0], object, mobile, target)

    except OSError:
        print ' Error: Cannot execute "%s", please provide full path to prosmart executable' % (exe)
        raise CmdException
    finally:
        shutil.rmtree(tempdir)

    if not quiet:
        print ' prosmart: done'

def _bfit_get_prior(distribution, em=0):
    from csb.statistics import scalemixture as sm

    if distribution == 'student':
        prior = sm.GammaPrior()
        if em: prior.estimator = sm.GammaPosteriorMAP()
    elif distribution == 'k':
        prior = sm.InvGammaPrior()
        if em: prior.estimator = sm.InvGammaPosteriorMAP()
    else:
        raise AttributeError('distribution')

    return prior

def xfit(mobile, target, mobile_state=-1, target_state=-1, load_b=0,
        cycles=10, match='align', guide=1, seed=0, quiet=1,
        bfit=0, distribution='student', _self=cmd):
    '''
DESCRIPTION

    Weighted superposition of the model in the first selection on to the model
    in the second selection. The weights are estimated with maximum likelihood.

    The result should be very similar to "theseus".

    Requires CSB, http://csb.codeplex.com

ARGUMENTS

    mobile = string: atom selection
 
    target = string: atom selection

    mobile_state = int: object state of mobile selection {default: current}

    target_state = int: object state of target selection {default: current}

    load_b = 0 or 1: save -log(weights) into B-factor column {default: 0}

SEE ALSO

    intra_xfit, align, super, fit, cealign, theseus
    '''
    from numpy import asarray, identity, log, dot, zeros
    from csb.bio.utils import distance_sq, wfit, fit
    from . import querying

    cycles, quiet = int(cycles), int(quiet)
    mobile_state, target_state = int(mobile_state), int(target_state)
    mobile_obj = querying.get_object_name(mobile, 1)

    if mobile_state < 1: mobile_state = querying.get_object_state(mobile_obj)
    if target_state < 1: target_state = querying.get_selection_state(target)

    if int(guide):
        mobile = '(%s) and guide' % (mobile)
        target = '(%s) and guide' % (target)

    mm = MatchMaker(mobile, target, match)

    Y = asarray(querying.get_coords(mm.mobile, mobile_state))
    X = asarray(querying.get_coords(mm.target, target_state))

    if int(seed):
        R, t = identity(3), zeros(3)
    else:
        R, t = fit(X, Y)

    if int(bfit):
        # adapted from csb.apps.bfit

        from csb.bio.utils import distance, probabilistic_fit
        from csb.statistics.scalemixture import ScaleMixture

        mixture = ScaleMixture(scales=X.shape[0],
                prior=_bfit_get_prior(distribution), d=3)

        for _ in range(cycles):
            data = distance(Y, dot(X - t, R))
            mixture.estimate(data)
            R, t = probabilistic_fit(X, Y, mixture.scales)

        scales = mixture.scales

    else:
        for _ in range(cycles):
            data = distance_sq(Y, dot(X - t, R))
            scales = 1.0 / data.clip(1e-3)
            R, t = wfit(X, Y, scales)

    m = identity(4)
    m[0:3,0:3] = R
    m[0:3,3] = t
    cmd.transform_object(mobile_obj, list(m.flat))

    if int(load_b):
        b_iter = iter(-log(scales))
        cmd.alter(mm.mobile, 'b = b_iter.next()', space=locals())

    if not quiet:
        print ' xfit: %d atoms aligned' % (len(X))

def intra_xfit(selection, load_b=0, cycles=20, guide=1, seed=0, quiet=1,
        bfit=0, distribution='student', _self=cmd):
    '''
DESCRIPTION

    Weighted superposition of all states of an object to the intermediate
    structure over all states. The weights are estimated with maximum
    likelihood.

    The result should be very similar to "intra_theseus".

    Requires CSB, http://csb.codeplex.com

ARGUMENTS

    selection = string: atom selection

    load_b = 0 or 1: save -log(weights) into B-factor column {default: 0}

NOTE

    Assumes all states to have identical number of CA-atoms.

SEE ALSO

    xfit, intra_fit, intra_theseus
    '''
    from numpy import asarray, identity, log, dot, zeros
    from csb.bio.utils import wfit, fit
    from .querying import get_ensemble_coords, get_object_name

    cycles, quiet = int(cycles), int(quiet)
    mobile_obj = get_object_name(selection, 1)
    n_models = cmd.count_states(mobile_obj)

    if int(guide):
        selection = '(%s) and guide' % (selection)

    X = asarray(get_ensemble_coords(selection))

    R, t = [identity(3)] * n_models, [zeros(3)] * n_models

    if int(bfit):
        # adapted from csb.apps.bfite

        from csb.bio.utils import average_structure, distance
        from csb.statistics.scalemixture import ScaleMixture

        average = average_structure(X)

        mixture = ScaleMixture(scales=X.shape[1],
                prior=_bfit_get_prior(distribution), d=3)

        for i in range(n_models):
            R[i], t[i] = fit(X[i], average)

        for _ in range(cycles):
            data = asarray([distance(average, dot(X[i] - t[i], R[i])) for i in range(n_models)])
            mixture.estimate(data.T)
            for i in range(n_models):
                R[i], t[i] = wfit(X[i], average, mixture.scales)

        scales = mixture.scales

    else:
        if int(seed):
            ensemble = X
        else:
            ensemble = []
            for i in range(n_models):
                R[i], t[i] = fit(X[i], X[0])
                ensemble.append(dot(X[i] - t[i], R[i]))

        for _ in range(cycles):
            ensemble = asarray(ensemble)
            average = ensemble.mean(0)
            data = ensemble.var(0).sum(1)
            scales = 1.0 / data.clip(1e-3)

            ensemble = []
            for i in range(n_models):
                R[i], t[i] = wfit(X[i], average, scales)
                ensemble.append(dot(X[i] - t[i], R[i]))

    m = identity(4)
    back = identity(4)
    back[0:3,0:3] = R[0]
    back[0:3,3] = t[0]

    for i in range(n_models):
        m[0:3,0:3] = R[i].T
        m[3,0:3] = -t[i]
        cmd.transform_object(mobile_obj, list(m.flat), state=i+1)

    # fit back to first state
    cmd.transform_object(mobile_obj, list(back.flat), state=0)

    if int(load_b):
        b_iter = iter(-log(scales))
        cmd.alter(selection, 'b = b_iter.next()', space=locals())

    if not quiet:
        print ' intra_xfit: %d atoms in %d states aligned' % (len(X[0]), n_models)

def promix(mobile, target, K=0, prefix=None, mobile_state=-1, target_state=-1,
        match='align', guide=1, quiet=1, async=-1, _self=cmd):
    '''
DESCRIPTION

    Finds rigid segments in two objects with different conformation.

    Requires CSB, http://csb.codeplex.com

ARGUMENTS

    mobile, target = string: atom selections

    K = integer: Number of segments {default: guess}

    prefix = string: Prefix of named segment selections to make

SEE ALSO

    intra_promix

REFERENCE

    Mixture models for protein structure ensembles
    Hirsch M, Habeck M. - Bioinformatics. 2008 Oct 1;24(19):2184-92
    '''
    from numpy import asarray
    from csb.statistics.mixtures import SegmentMixture as Mixture
    from .querying import get_coords, get_object_name

    K, guide, quiet, async = int(K), int(guide), int(quiet), int(async)
    mobile_state, target_state = int(mobile_state), int(target_state)
    if async < 0:
        async = not quiet

    if isinstance(target, str) and target.isdigit() and \
            cmd.count_atoms('?' + target) == 0 and cmd.count_states(mobile) > 1:
        print ' Warning: sanity test suggest you want "intra_promix"'
        return intra_promix(mobile, target, prefix, 0, guide, quiet, async)

    if guide:
        mobile = '(%s) and guide' % (mobile)
        target = '(%s) and guide' % (target)

    cmd.color('gray', mobile)
    obj = get_object_name(mobile)
    mm = MatchMaker(mobile, target, match)
    selection = mm.mobile

    X = asarray([
        get_coords(mm.mobile, mobile_state),
        get_coords(mm.target, target_state),
    ])

    if not async:
        _promix(**locals())
    else:
        import threading
        t = threading.Thread(target=_promix, kwargs=locals())
        t.setDaemon(1)
        t.start()

def intra_promix(selection, K=0, prefix=None, conformers=0, guide=1,
        quiet=1, async=-1, _self=cmd):
    '''
DESCRIPTION

    Finds rigid segments in a multi-state object.

    Requires CSB, http://csb.codeplex.com

ARGUMENTS

    selection = string: atom selection

    K = integer: Number of segments {default: guess}

    prefix = string: Prefix of named segment selections to make

SEE ALSO

    promix

REFERENCE

    Mixture models for protein structure ensembles
    Hirsch M, Habeck M. - Bioinformatics. 2008 Oct 1;24(19):2184-92
    '''
    from numpy import asarray
    from csb.statistics import mixtures
    from .querying import get_ensemble_coords, get_object_name

    K, conformers = int(K), int(conformers)
    guide, quiet, async = int(guide), int(quiet), int(async)
    if async < 0:
        async = not quiet

    Mixture = mixtures.ConformerMixture if conformers else mixtures.SegmentMixture

    obj = get_object_name(selection)
    n_models = cmd.count_states(obj)

    if guide:
        selection = '(%s) and guide' % (selection)

    if n_models < 2:
        print ' Error: object needs multiple states'
        raise CmdException

    X = asarray(get_ensemble_coords(selection))
    assert X.shape == (n_models, cmd.count_atoms(selection), 3)

    if not async:
        _promix(**locals())
    else:
        import threading
        t = threading.Thread(target=_promix, kwargs=locals())
        t.setDaemon(1)
        t.start()

def _promix(**kwargs):
    conformers = 0

    exec ','.join(kwargs) + ', = kwargs.values()'

    if not prefix:
        if conformers:
            prefix = obj + '_conformer'
        else:
            prefix = obj + '_segment'
    cmd.delete(prefix + '_*')

    id_list = []
    cmd.iterate(selection, 'id_list.append(ID)', space=locals())

    mixture = Mixture.new(X, K)
    membership = mixture.membership

    if conformers:
        states_list = [0] * mixture.K
        for (i,k) in enumerate(membership):
            states_list[k] += 1
            name = '%s_%d' % (prefix, k+1)
            cmd.create(name, obj, i+1, states_list[k])
    else:
        cmd.color('gray', selection)
        for k in range(mixture.K):
            name = '%s_%d' % (prefix, k+1)
            id_list_k = [i for (i, m) in zip(id_list, membership) if m == k]
            cmd.select_list(name, obj, id_list_k)
            cmd.disable(name)
            cmd.color(k + 2, name)

    for k, (sigma, w) in enumerate(zip(mixture.sigma, mixture.w)):
        print ' %s_%d: sigma = %6.3f, w = %.3f' % (prefix, k+1, sigma, w)

    print ' BIC: %.2f' % (mixture.BIC)
    print ' Log Likelihood: %.2f' % (mixture.log_likelihood)

# all those have kwargs: mobile, target, mobile_state, target_state
align_methods = ['align', 'super', 'cealign', 'tmalign', 'theseus', 'prosmart', 'xfit']
align_methods_sc = cmd.Shortcut(align_methods)

# pymol commands
cmd.extend('alignwithanymethod', alignwithanymethod)
cmd.extend('tmalign', tmalign)
cmd.extend('dyndom', dyndom)
cmd.extend('gdt_ts', gdt_ts)
cmd.extend('local_rms', local_rms)
cmd.extend('extra_fit', extra_fit)
cmd.extend('intra_theseus', intra_theseus)
cmd.extend('theseus', theseus)
cmd.extend('prosmart', prosmart)
cmd.extend('xfit', xfit)
cmd.extend('intra_xfit', intra_xfit)
cmd.extend('promix', promix)
cmd.extend('intra_promix', intra_promix)

# autocompletion
_auto_arg0_align = cmd.auto_arg[0]['align']
_auto_arg1_align = cmd.auto_arg[1]['align']
cmd.auto_arg[0].update([
    ('alignwithanymethod', _auto_arg0_align),
    ('tmalign', _auto_arg0_align),
    ('dyndom', _auto_arg0_align),
    ('gdt_ts', _auto_arg0_align),
    ('local_rms', _auto_arg0_align),
    ('extra_fit', _auto_arg0_align),
    ('theseus', _auto_arg0_align),
    ('intra_theseus', _auto_arg1_align),
    ('prosmart', _auto_arg0_align),
    ('xfit', _auto_arg0_align),
    ('intra_xfit', _auto_arg0_align),
    ('promix', _auto_arg0_align),
    ('intra_promix', _auto_arg0_align),
])
cmd.auto_arg[1].update([
    ('alignwithanymethod', _auto_arg1_align),
    ('tmalign', _auto_arg1_align),
    ('dyndom', _auto_arg1_align),
    ('gdt_ts', _auto_arg1_align),
    ('local_rms', _auto_arg1_align),
    ('extra_fit', cmd.auto_arg[0]['disable']),
    ('theseus', _auto_arg1_align),
    ('prosmart', _auto_arg1_align),
    ('xfit', _auto_arg1_align),
    ('promix', _auto_arg0_align),
])
cmd.auto_arg[2].update([
    ('extra_fit', [ align_methods_sc, 'alignment method', '' ]),
])

# vi: ts=4:sw=4:smarttab:expandtab
