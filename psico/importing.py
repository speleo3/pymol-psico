'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

import os.path
from pymol import cmd, CmdException

mysql_kwargs = {
    'host': 'abt1-th.eb.local',
    'user': 'bioguest',
    'passwd': '',
    'db': 'biomol',
}

def local_mirror_pdb(code):
    '''
DESCRIPTION

    API only. Returns the filename for given PDB code. If you have a local PDB
    mirror, you probably want to overwrite this function.
    '''
    return '/ebio/abt1_toolkit/share/wye/databases/pdb/all/pdb%s.ent' % (code.lower())

def fetch(code, name='', type='pdb', quiet=1, async=-1, **kwargs):
    '''
PSICO NOTES

    This is an enhanced version of the "fetch" command that also handles
    SCOP and CATH identifiers.

    You can also fetch from a local filesystem mirror by setting up a
    function "psico.importing.local_mirror_pdb(code)" that returns a
    full path to a pdb file. 
    '''
    from pymol.importing import fetch as pymol_fetch

    quiet, async = int(quiet), int(async)

    if type != 'pdb':
        return pymol_fetch(code, name, type=type, quiet=quiet, async=async, **kwargs)

    load_kwargs = {'state':0, 'format':'', 'finish':1, 'discrete':-1,
            'quiet':quiet, 'multiplex':None, 'zoom':-1, 'partial':0, 'mimic':1}
    load_kwargs = dict((k,kwargs[k]) for k in kwargs if k in load_kwargs)
    code_list = code.split()
    for code in list(code_list):
        if len(code) in [5,6]:
            fetch_chain(code, name, **kwargs)
            code_list.remove(code)
        elif len(code) == 7:
            if code[0] == 'd':
                fetch_scop(code, name, **kwargs)
            elif code[0].isdigit():
                fetch_cath(code, name, **kwargs)
            else:
                print ' Warning: 7 character ID but does not look like SCOP or CATH'
                print ' Warning: skipping', code
            code_list.remove(code)
        else:
            filename = local_mirror_pdb(code)
            if os.path.exists(filename):
                if not quiet:
                    print ' Using local mirror:', filename
                cmd.load(filename, name if len(name) else code, **load_kwargs)
                code_list.remove(code)
    if len(code_list) > 0:
        return pymol_fetch(' '.join(code_list), name, quiet=quiet, async=async, **kwargs)
fetch.__doc__ += cmd.fetch.__doc__

cath_domains = {}
def cath_parse_domall(filename=''):
    ''' 
    Download "CathDomall.seqreschopping" to fetch_path (once), parse it and
    store results in global dict "cath_domains".
    '''
    if filename == '':
        fetch_path = cmd.get('fetch_path')
        if fetch_path == '.':
            import tempfile
            fetch_path = tempfile.gettempdir()
        basename = 'CathDomall.seqreschopping'
        filename = os.path.join(fetch_path, basename)
        if not os.path.exists(filename):
            import urllib2
            handle = urllib2.urlopen('http://release.cathdb.info/latest/' + basename)
            out = open(filename, 'w')
            out.write(handle.read())
            out.close()

    for line in open(os.path.expanduser(filename)):
        domid, resix = line.split()
        domain = []
        for resiy in resix.split(','):
            resi = resiy.split('-')
            domain.append((domid[4], resi[0], resi[1]))
        cath_domains[domid] = domain

def fetch_cath(code, name='', **kwargs):
    if name == '':
        name = code
    r = fetch(code[:4], name, async=0, **kwargs)
    if cmd.is_error(r):
        return r
    if code[5:7] == '00':
        if code[4] != '0':
            cmd.remove(name + ' and not chain ' + code[4])
    else:
        try:
            if len(cath_domains) == 0:
                cath_parse_domall()
            d = cath_domains[code]
            rsele = ''
            for frag in d:
                chain, start, stop = frag
                rsele += ' (chain ' + chain
                if start and stop:
                    rsele += ' and resi %s-%s' % (start, stop)
                rsele += ')'
            cmd.remove(name + ' and not (' + rsele + ')')
        except:
            print ' Warning: CATH domain resiude range handling failed'
            return -1

def fetch_scop(code, name='', **kwargs):
    if name == '':
        name = code
    r = fetch(code[1:5], name, async=0, **kwargs)
    if cmd.is_error(r):
        return r
    try:
        import Bio.SCOP
        import MySQLdb
        conn = MySQLdb.connect(**mysql_kwargs)
        scop = Bio.SCOP.Scop(db_handle=conn)
        d = scop.getDomainBySid(code)
        if len(d.residues.fragments) > 0:
            rsele = ''
            for frag in d.residues.fragments:
                chain, start, stop = frag
                rsele += ' (chain ' + chain
                if start and stop:
                    rsele += ' and resi %s-%s' % (start, stop)
                rsele += ')'
            cmd.remove(name + ' and not (' + rsele + ')')
    except:
        print ' Warning: SCOP domain resiude range handling failed'
        return -1

def fetch_chain(code, name='', **kwargs):
    if name == '':
        name = code
    chain = code[4] if len(code) == 5 else code[5]
    r = fetch(code[:4], name, async=0, **kwargs)
    if cmd.is_error(r):
        return r
    cmd.remove(name + ' and not chain ' + chain)

def loadall(pattern, group='', quiet=1, **kwargs):
    '''
DESCRIPTION

    Load all files matching given globbing pattern
    '''
    import glob, os
    filenames = glob.glob(cmd.exp_path(pattern))
    for filename in filenames:
        if not quiet:
            print ' Loading', filename
        cmd.load(filename, **kwargs)
    if len(group):
        if kwargs.get('object', '') != '':
            print ' Warning: group and object arguments given'
            members = [kwargs['object']]
        else:
            from pymol.cmd import file_ext_re, gz_ext_re, safe_oname_re
            members = [gz_ext_re.sub('', file_ext_re.sub('',
                safe_oname_re.sub('_', os.path.split(filename)[-1])))
                for filename in filenames]
        cmd.group(group, ' '.join(members))

def load_traj_dcd(filename, object='', state=0, start=1, stop=-1, quiet=1):
    '''
DESCRIPTION

    Load DCD trajectory.

    http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

SEE ALSO

    save_traj, load_traj_crd, load_traj
    '''
    import struct
    import itertools

    state, quiet = int(state), int(quiet)
    start, stop = int(start), int(stop)

    if object == '':
        object = os.path.splitext(os.path.basename(filename))[0]
    if object not in cmd.get_object_list():
        print ' Error: must load object topology before loading trajectory'
        raise CmdException

    if state < 1:
        state = cmd.count_states(object) + 1
    natom = cmd.count_atoms(object)

    endian = '<'
    handle = open(filename)
    read = lambda fmt: struct.unpack(endian + fmt,
            handle.read(struct.calcsize(fmt)))

    # Header
    length = read('i')[0]
    if length != 0x54:
        if not quiet:
            print ' Info: big-endian'
        endian = '>'

    fmt = '4s i i i 5i i d 9i i'
    HDR, NSET, ISTRT, NSAVC, _, _, _, _, _, NATOMNFREAT, DELTA, \
            _, _, _, _, _, _, _, _, _, LENGTH = read(fmt)
    assert HDR == 'CORD'
    assert LENGTH == length

    # Title
    length = read('i')[0]
    fmt = 'i %ds i' % (length - 4)
    NTITLE, TITLE, LENGTH = read(fmt)
    assert LENGTH == length

    # NATOM
    length, NATOM, LENGTH = read('iii')
    assert LENGTH == length
    assert NATOM == natom

    # Coord Sets
    fmt = 'i %df i' % (NATOM)
    def readX():
        X = read(fmt)
        assert X[0] == X[-1]
        return X[1:-1]

    for frame in itertools.count(1):
        if 0 < stop < frame or 0 < NSET < frame:
            break

        length = read('i')[0]
        if length == 6*8:
            box = read('6d i')
            assert box[-1] == length
        else:
            handle.seek(-struct.calcsize('i'), 1)

        try:
            XYZ = readX(), readX(), readX()
        except struct.error:
            break

        if frame < start:
            continue

        coordset = map(list, zip(*XYZ))
        assert len(coordset) == natom

        load_coords(coordset, object, state)
        state += 1

def load_traj_crd(filename, object='', state=0, box=0, start=1, stop=-1, quiet=1):
    '''
DESCRIPTION

    Load Amber CRD trajectory file.

    http://ambermd.org/formats.html#trajectory

    For some (unkown to me) reason, this does not work with load_traj (but it
    should!).

SEE ALSO

    load_traj
    '''
    state, box, quiet = int(state), int(box), int(quiet)
    start, stop = int(start), int(stop)
    if object == '':
        object = os.path.splitext(os.path.basename(filename))[0]
    if object not in cmd.get_object_list():
        print ' Error: must load object topology before loading trajectory'
        raise CmdException
    if state < 1:
        state = cmd.count_states(object) + 1
    natom = cmd.count_atoms(object)
    line_it = iter(open(filename))
    line_it.next() # skip title
    def crd_coord_iter():
        '''
        Iterator that yields True at the beginning of a coord set, followed
        by the sequence of coordinates (floats).

        6F12.7 (rst7)
        10F8.3 (crd)
        '''
        count = 0
        ncoord = natom * 3
        width = 0
        for line in line_it:
            if width == 0:
                # in second line look for floating number width
                second_dot = line.find('.', 8)
                if second_dot == -1:
                    if not quiet:
                        print 'Number of atoms in second line'
                    # assume number of atoms in second line
                    _natom = int(line)
                    assert _natom == natom, 'Numbers differ: %d, %d' % (natom, _natom)
                    line = line_it.next()
                    second_dot = line.find('.', 8)
                if second_dot == 12:
                    width = 8
                elif second_dot == 16:
                    width = 12
                else:
                    print 'WARNING: could neither detect 6F12.7 nor 10F8.3 format'
                    raise CmdException
            if count > ncoord:
                print 'Error!', ncoord, count
                raise CmdException
            if count == ncoord:
                count = 0
                if box:
                    if not quiet:
                        print 'skipping box'
                    # skip box line
                    continue
            for i in range(0, len(line.rstrip()), width):
                x = float(line[i:i+width])
                if count == 0:
                    yield True
                count += 1
                yield x
        yield False
    frame = 1
    coord_it = crd_coord_iter()
    while coord_it.next() == True:
        coordset = [[coord_it.next(), coord_it.next(), coord_it.next()]
                for _ in range(natom)]
        if frame >= start:
            load_coords(coordset, object, state)
            state += 1
        if stop > 0 and frame == stop:
            break
        frame += 1

def load_3d(filename, object=''):
    '''
DESCRIPTION

    Load a survex 3d cave survey as "molecule"

    http://survex.com
    http://trac.survex.com/browser/trunk/doc/3dformat.htm
    '''
    from chempy import Atom, Bond, models
    from struct import unpack

    if object == '':
        object = os.path.splitext(os.path.basename(filename))[0]

    f = open(filename, 'rb')

    line = f.readline() # File ID
    if not line.startswith('Survex 3D Image File'):
        print " Error: not a Survex 3D File"
        raise CmdException

    line = f.readline() # File format version
    assert line[0] == 'v'
    ff_version = int(line[1:])

    line = unicode(f.readline(), 'latin1') # Survex title
    line = f.readline() # Timestamp

    class Station:
        def __init__(self):
            self.labels = []
            self.adjacent = []
            self.lrud = None
            self.flag = 0
        def connect(self, other):
            self.adjacent.append(other)
        def is_surface(self):
            return self.flag & 0x01
        def is_underground(self):
            return self.flag & 0x02
        def is_entrance(self):
            return self.flag & 0x04
        def is_exported(self):
            return self.flag & 0x08
        def is_fixed(self):
            return self.flag & 0x10

    class Survey(dict):
        def __init__(self):
            self.prev = None
            self.curr_label = ''
            self.labelmap = {}
        def get(self, xyz):
            return dict.setdefault(self, tuple(xyz), Station())
        def line(self, xyz):
            s = self.get(xyz)
            self.prev.connect(s)
            self.prev = s
        def move(self, xyz):
            s = self.get(xyz)
            self.prev = s
        def label(self, xyz, flag=0):
            s = self.get(xyz)
            s.labels.append(self.curr_label)
            self.labelmap[s.labels[-1]] = s
            if flag > 0:
                s.flag = flag
        def lrud(self, lrud):
            s = self.labelmap[self.curr_label]
            s.lrud = lrud

    survey = Survey()

    def read_xyz():
        return unpack('<iii', f.read(12))

    def read_len():
        len = read_byte()
        if len == 0xfe:
            len += unpack('<H', f.read(2))[0]
        elif len == 0xff:
            len += unpack('<I', f.read(4))[0]
        return len

    def read_label():
        len = read_len()
        if len > 0:
            survey.curr_label += skip_bytes(len)

    def skip_bytes(n):
        return f.read(n)

    def read_byte():
        byte = f.read(1)
        if len(byte) != 1:
            return -1
        return ord(byte)

    while 1:
        byte = read_byte()
        if byte == -1:
            break
        
        if byte == 0x00:
            # STOP
            survey.curr_label = ''
        elif byte <= 0x0e:
            # TRIM
            # FIXME: according to doc, trim 16 bytes, but img.c does 17!
            (i,n) = (-17,0)
            while n < byte:
                i -= 1
                if survey.curr_label[i] == '.': n += 1
            survey.curr_label = survey.curr_label[:i + 1]
        elif byte <= 0x0f:
            # MOVE
            xyz = read_xyz()
            survey.move(xyz)
        elif byte <= 0x1f:
            # TRIM
            survey.curr_label = survey.curr_label[:15 - byte]
        elif byte <= 0x20:
            # DATE
            if ff_version < 7:
                skip_bytes(4)
            else:
                skip_bytes(2)
        elif byte <= 0x21:
            # DATE
            if ff_version < 7:
                skip_bytes(8)
            else:
                skip_bytes(3)
        elif byte <= 0x22:
            # Error info
            skip_bytes(5 * 4)
        elif byte <= 0x23:
            # DATE
            skip_bytes(4)
        elif byte <= 0x24:
            # DATE
            continue
        elif byte <= 0x2f:
            # Reserved
            continue
        elif byte <= 0x31:
            # XSECT
            read_label()
            lrud = unpack('<hhhh', f.read(8))
            survey.lrud(lrud)
        elif byte <= 0x33:
            # XSECT
            read_label()
            lrud = unpack('<iiii', f.read(16))
            survey.lrud(lrud)
        elif byte <= 0x3f:
            # Reserved
            continue
        elif byte <= 0x7f:
            # LABEL
            read_label()
            xyz = read_xyz()
            survey.label(xyz, byte & 0x3f)
        elif byte <= 0xbf:
            # LINE
            read_label()
            xyz = read_xyz()
            survey.line(xyz)
        elif byte <= 0xff:
            # Reserved
            continue

    model = models.Indexed()
    for (xyz,s) in survey.iteritems():
        l0, _, l1 = s.labels[0].rpartition('.')
        resi, name = l1[:5], l1[5:]
        segi, chain, resn = l0[-8:-4], l0[-4:-3], l0[-3:]
        atom = Atom()
        atom.coord = [i/100.0 for i in xyz]
        atom.segi = segi
        atom.chain = chain
        atom.resn = resn
        atom.name = name
        atom.resi = resi
        atom.b = atom.coord[2]
        atom.label = s.labels[0]
        if s.lrud is not None:
            atom.vdw = sum(s.lrud)/400.0
        model.add_atom(atom)

    s2i = dict((s,i) for (i,s) in enumerate(survey.itervalues()))
    for (s,i) in s2i.iteritems():
        for o in s.adjacent:
            bnd = Bond()
            bnd.index = [i, s2i[o]]
            model.add_bond(bnd)

    cmd.load_model(model, object, 1)
    cmd.show_as('lines', object)
    cmd.spectrum('b', 'rainbow', object)

def set_raw_alignment(name, aln, transform=0):
    '''
DESCRIPTION

    API only.
    Load an alignment object from a list like the one obtained with
    cmd.get_raw_alignment
    '''
    import itertools
    if not isinstance(aln[0], dict):
        aln = [dict(idx_pair) for idx_pair in aln]
    models = set(model for idx_pdict in aln for model in idx_pdict)
    sele1 = cmd.get_unused_name('_sele1')
    sele2 = cmd.get_unused_name('_sele2')
    fit = cmd.fit if transform else cmd.rms_cur
    for model1, model2 in itertools.combinations(models, 2):
        index_list1 = []
        index_list2 = []
        for idx_pdict in aln:
            if model1 in idx_pdict and model2 in idx_pdict:
                index_list1.append(idx_pdict[model1])
                index_list2.append(idx_pdict[model2])
        cmd.select_list(sele1, model1, index_list1, mode='index')
        cmd.select_list(sele2, model2, index_list2, mode='index')
        fit(sele1, sele2, cycles=0, matchmaker=4, object=name)
    cmd.delete(sele1)
    cmd.delete(sele2)

def load_aln(filename, object=None, mobile=None, target=None, mobile_id=None,
        target_id=None, format='', transform=0, quiet=1):
    '''
DESCRIPTION

    Load an alignment from file and apply it to already loaded structures.

USAGE

    load_aln filename [, object [, mobile [, target [, mobile_id [, target_id [, format ]]]]]]

ARGUMENTS

    filename = string: alignment file

    object = string: name of the object {default: filename prefix}

    mobile, target = string: atom selections {default: ids from alignment file}

    mobile_id, target_id = string: ids from alignment file {default: first two}

    format = string: file format, see http://biopython.org/wiki/AlignIO
    {default: guess from first line in file}

EXAMPLE

    fetch 1bz4 1cpr, async=0
    super 1bz4 and guide, 1cpr and guide, object=aln1, window=5
    save /tmp/super.aln, aln1
    delete aln1
    load_aln /tmp/super.aln, aln2, format=clustal
    '''
    from . import one_letter
    from .seqalign import needle_alignment, alignment_mapping, aln_magic_read

    quiet = int(quiet)
    if object is None:
        object = os.path.basename(filename).rsplit('.', 1)[0]

    # load alignment file
    alignment = aln_magic_read(cmd.exp_path(filename), format)
    aln_dict = dict((record.id, record) for record in alignment)
    mobile_record = alignment[0] if mobile_id is None else aln_dict[mobile_id]
    target_record = alignment[1] if target_id is None else aln_dict[target_id]

    # guess selections from sequence identifiers (if not given)
    if mobile is None: mobile = mobile_record.id
    if target is None: target = target_record.id

    try:
        mobile_obj = cmd.get_object_list('(' + mobile + ')')[0]
        target_obj = cmd.get_object_list('(' + target + ')')[0]
    except:
        print ' Error: selection "%s" or "%s" does not exist' % (mobile, target)
        raise CmdException

    # get structure models and sequences
    mobile_model = cmd.get_model('(%s) and guide' % mobile)
    target_model = cmd.get_model('(%s) and guide' % target)
    mobile_sequence = ''.join(one_letter.get(a.resn, 'X') for a in mobile_model.atom)
    target_sequence = ''.join(one_letter.get(a.resn, 'X') for a in target_model.atom)

    # align sequences from file to structures
    mobile_aln = needle_alignment(str(mobile_record.seq), mobile_sequence)
    target_aln = needle_alignment(str(target_record.seq), target_sequence)

    # get index mappings
    mobile_aln_i2j = dict(alignment_mapping(*mobile_aln))
    target_aln_i2j = dict(alignment_mapping(*target_aln))
    record_i2j = alignment_mapping(mobile_record, target_record)

    # build alignment list
    r = []
    for i, j, in record_i2j:
        if i in mobile_aln_i2j and j in target_aln_i2j:
            i = mobile_aln_i2j[i]
            j = target_aln_i2j[j]
            r.append({
                mobile_obj: mobile_model.atom[i].index,
                target_obj: target_model.atom[j].index,
            })

    set_raw_alignment(object, r, int(transform))
    return r

def load_gro(filename, object='', state=-1, quiet=1, zoom=-1):
    '''
DESCRIPTION

    Load a gro file (molecular structure in Gromos87 format).

    http://manual.gromacs.org/online/gro.html

    This command can handle larger ID and resi values than PDB files and thus
    is better (but slower) than converting a gro file to PDB format using
    editconf.

    Based on a script from Tsjerk Wassenaar, posted on pymol-users mailing list.
    http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg09394.html
    '''
    if object == '':
        object = os.path.basename(filename).rsplit('.', 1)[0]

    _pdbAtomLine = 'ATOM  %5i  %-3s %3s%2s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s  '
    _pdbBoxLine = 'CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1'

    def pdbAtom(line):
        ''' Parses a gro atom line and returns a pdb ATOM line '''
        resi = int(line[:5])
        resn = line[5:10].strip()
        name = line[10:15].strip()
        ID = int(line[15:20])
        x, y, z = 10*float(line[20:28]), 10*float(line[28:36]), 10*float(line[36:44])
        segi = '%02d%02d' % (ID/100000, resi/10000)
        ID, resi = ID%100000, resi%10000
        return _pdbAtomLine % (ID, name[:3], resn[:3], 'A', resi, x, y, z, 1, 0, segi, name[0])

    def pdbBox(line):
        ''' Parses a gro box line and returns a pdb CRYST1 line '''
        from math import degrees
        from chempy.cpv import length, get_angle

        v = [10*float(i) for i in line.split()] + 6*[0] # Padding for rectangular boxes
        v1, v2, v3 = (v[0], v[3], v[4]), (v[5], v[1], v[6]), (v[7], v[8], v[2])

        a = length(v1)
        b = length(v2)
        c = length(v3)

        alpha = degrees(get_angle(v2, v3))
        beta = degrees(get_angle(v1, v3))
        gamma = degrees(get_angle(v1, v2))

        return _pdbBoxLine % (a, b, c, alpha, beta, gamma)

    pdb = ['']
    stream = open(filename)
    for model in xrange(1, 9999):
        title = stream.readline()
        natoms = stream.readline().strip()
        if natoms == '':
            break
        natoms = int(natoms)

        pdb.append('MODEL %8d' % model)
        for i in range(natoms):
            pdb.append(pdbAtom(stream.readline()))
        pdb.append('ENDMDL')

        boxline = stream.readline()
        if model == 1:
            pdb[0] = pdbBox(boxline)
    stream.close()
    pdb.extend(['END', ''])

    cmd.read_pdbstr('\n'.join(pdb), object, int(state), quiet=int(quiet), zoom=int(zoom))
    cmd.alter(object, '(ID,resi,segi) = (ID+100000*int(segi[:2]),resv+10000*int(segi[2:]),"")')

def load_coords(coords, object, state, quiet=1):
    '''
DESCRIPTION

    API only. Load object coordinates.

    TODO: How to handle the object matrix?

SEE ALSO

    cmd.load_coordset (new in 1.7.3.0)
    pymol.experimenting.load_coords (considered broken!)
    '''
    if hasattr(cmd, 'load_coordset'):
        return cmd.load_coordset(coords, object, int(state))

    if not (isinstance(coords, list) and isinstance(coords[0], list)):
        coords = map(list, coords)
    r = cmd._cmd.load_coords(cmd._COb, object, coords, int(state)-1,
            cmd.loadable.model)
    if cmd._raising(r): raise CmdException
    return r

def load_mtz(filename, prefix='', maptypes='FoFc 2FoFc', multistate=0, quiet=1):
    '''
DESCRIPTION

    Load a MTZ file as two map objects (FoFc, 2FoFc)

    This only works with the incentive PyMOL product!

USAGE

    load_mtz filename [, prefix [, maptypes ]]
    '''
    from pymol import headering

    multistate, quiet = int(multistate), int(quiet)

    header = headering.MTZHeader(filename)

    for coltype in ('F', 'P'):
        if not header.getColumnsOfType(coltype):
            print ' Error: could not find %s type column' % coltype
            raise CmdException

    if not prefix:
        prefix = os.path.basename(filename).rsplit('.', 1)[0]

    for maptype in maptypes.split():
        F, P, prg = header.guessCols(maptype)
        if None in (F, P):
            print ' Warning: No %s map' % (maptype)
            continue

        name = prefix + '.' + maptype
        if not multistate:
            cmd.delete(name)

        cmd.map_generate(name, filename, F, P, 'None', 0, 0, quiet)
        if name not in cmd.get_names('objects'):
            print ' Error: Loading %s map failed.' % (maptype)
            print ' This PyMOL version might not be capable of loading MTZ files'
            raise CmdException

# commands
cmd.extend('loadall', loadall)
cmd.extend('load_traj_crd', load_traj_crd)
cmd.extend('load_traj_dcd', load_traj_dcd)
cmd.extend('load_3d', load_3d)
cmd.extend('load_aln', load_aln)
cmd.extend('load_gro', load_gro)
cmd.extend('load_mtz', load_mtz)

# vi:expandtab:smarttab
