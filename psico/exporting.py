'''
File export module that provides overloaded save commands with secondary
structure and crystal information header, as well as saving to trajectory
formats.

(c) 2010-2012 Thomas Holder and Steffen Schmidt, MPI for Developmental Biology
(c) 2009 Sean Law, Michigan State University (save2traj)

License: BSD-2-Clause
'''

from __future__ import print_function

try:
    file
except NameError:
    from io import FileIO as file

from pymol import cmd, CmdException
from pymol import selector

def _assert_package_import():
    if not __name__.endswith('.exporting'):
        raise CmdException("Must do 'import psico.exporting' instead of 'run ...'")

## trajectory stuff

def save_traj(filename, selection='(all)', format='', box=0, quiet=1):
    '''
DESCRIPTION

    Save coordinates of a multi-state object to a trajectory file (DCD OR CRD).

    Based on http://pymolwiki.org/index.php/Save2traj by Sean Law

USAGE

    save_traj filename [, selection [, format ]]

ARGUMENTS

    filename = string: file path to be written

    selection = string: atoms to save

    format = string: 'dcd' or 'crd' (alias 'charmm' or 'amber') {default:
    determined from filename extension)
    '''
    _assert_package_import()
    from . import querying

    box = int(box)

    # Get NATOMS, NSTATES
    NATOMS = cmd.count_atoms(selection)
    NSTATES = cmd.count_states(selection)

    # Determine Trajectory Format
    if format == '' and '.' in filename:
        format = filename.rsplit('.', 1)[1]
    format = format.lower()
    if format in ['charmm', 'dcd']:
        Outfile = DCDOutfile
    elif format in ['amber', 'trj', 'crd']:
        Outfile = CRDOutfile
    elif format in ['rst', 'rst7']:
        Outfile = RSTOutfile
    else:
        print('Unknown format:', format)
        raise CmdException

    f = open(filename, 'wb')

    # size of periodic box
    if box:
        try:
            boxdim = cmd.get_symmetry(selection)[0:3]
        except:
            boxdim = [0,0,0]
    else:
        boxdim = None

    outfile = Outfile(filename, NSTATES, NATOMS, box=boxdim)

    # Write Trajectory Coordinates
    for state in range(1, NSTATES+1):
        xyz = querying.get_coords(selection, state)
        outfile.writeCoordSet(xyz)

    outfile.close()

    if not int(quiet):
        fmt = 'Wrote trajectory in %s format with %d atoms and %d frames to file %s'
        print(fmt % (format, NATOMS, NSTATES, filename))

class DCDOutfile(file):
    '''
http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
    '''
    def __init__(self, filename, nstates, natoms, vendor='PyMOL', box=None):
        file.__init__(self, filename, 'wb')
        self.natoms = natoms
        self.fmt = '%df' % (natoms)
        charmm = int(nstates > 0)

        # Header
        fmt='4s 9i d 9i'
        header = [b'CORD', # 4s
                nstates, 1, 1, 0, 0, 0, 0, 0, 0, # 9i
                1.0, # d
                0, 0, 0, 0, 0, 0, 0, 0, 0, # 9i
                ]
        if charmm:
            # DELTA is stored as a double with X-PLOR but as a float with CHARMm
            fmt = '4s 9i f 10i'
            header.append(24) # dummy charmm version number
        self.writeFortran(header,fmt)
     
        # Title
        fmt = 'i80s80s'
        title = [2, # 1i
                b'* TITLE'.ljust(80), # 80s
                (b'* Created by ' + vendor.encode()).ljust(80), # 80s
                ]
        self.writeFortran(title,fmt,length=160+4)
 
        # NATOM
        self.writeFortran([natoms],'i')

    def writeFortran(self, buffer, fmt, length=0):
        '''
        Write FORTRAN unformatted binary record.
        '''
        import struct
        if length == 0:
            length = struct.calcsize(fmt)
        self.write(struct.pack('i', length))
        self.write(struct.pack(fmt, *buffer))
        self.write(struct.pack('i', length))

    def writeCoordSet(self, xyz, transposed=0):
        '''
        Write a 3xNATOMS coord matrix.
        '''
        if not transposed:
            xyz = list(zip(*xyz))
        assert len(xyz) == 3, 'Wrong number of dimensions'
        for coor in xyz:
            assert len(coor) == self.natoms, 'Wrong number of atoms'
            self.writeFortran(coor, self.fmt)

class CRDOutfile(file):
    '''
http://ambermd.org/formats.html#trajectory
    '''
    fmt = '%8.3f'
    columns = 10

    def __init__(self, filename, nstates=-1, natoms=-1, vendor='PyMOL', box=None):
        file.__init__(self, filename, 'w')
        self.natoms = natoms
        self.box = box

        # Write Trajectory Header Information
        print('TITLE : Created by %s with %d atoms' % (vendor, natoms), file=self)

    def writeCoordSet(self, xyz, transposed=0):
        '''
        Write a NATOMSx3 coord matrix.
        '''
        if transposed:
            xyz = list(zip(*xyz))
        if self.natoms > -1:
            assert len(xyz) == self.natoms, 'Wrong number of atoms'
        assert len(xyz[0]) == 3, 'Wrong number of dimensions'
        f = self
        count = 0
        for coord in xyz:
            for c in coord:
                f.write(self.fmt % c)
                count += 1
                if count % self.columns == 0:
                    f.write('\n')
        if count % self.columns != 0:
            f.write('\n')

        # size of periodic box
        if self.box is not None:
            for c in self.box:
                f.write(self.fmt % c)
            f.write('\n')

class RSTOutfile(CRDOutfile):
    '''
http://ambermd.org/formats.html#restart
    '''
    fmt = '%12.7f'
    columns = 6

    def __init__(self, *args, **kwargs):
        super(RSTOutfile, self).__init__(*args, **kwargs)

        print('%5i%s' % (self.natoms, '  0.0000000e+00' * 5), file=self)

## pdb header stuff

def get_pdb_sss(selection='(all)', state=-1, quiet=1):
    '''
DESCRIPTION

    API-only. Return the PDB "Secondary Structure Section" for a given
    selection to put in the header section of a PDB file. Takes the "ss"
    atom property of CA atoms.

    http://www.wwpdb.org/documentation/format33/sect5.html

ARGUMENT

    selection = string: atom selection

    state = int: object state {default: -1}
    '''
    ss = {}         # storing the secondary structure elements by
                    # {chain}{ss.type}[ start_atom_object,
                    # end_atom_object ]

    # Get a list of CA atoms and read the secondary structure
    # annotation This loop assumes that the atoms are in consecutive
    # order i.e. sorted by chain & resi
    for at in cmd.get_model( '(' + selection + ') and n. CA and polymer',
                             state=state).atom:
        if at.ss == '': continue 

        # Init ss dictionary if key / ss doesn't exist
        L = ss.setdefault(at.chain, {}).setdefault(at.ss, [])

        # Check if a new ss has to be expanded (replace last atom in ss with at)
        # or else a new ss will be appended
        if len(L) and L[-1][1].resi_number == (at.resi_number - 1):
            L[-1][1] = at
        else:
            L.append([at, at])

    ssstr = []          # the output string

    # Iterate over stored secondary structures and add formatted
    # string to ssstr
    for chain in ss:
        for s in ss[chain]:
            for i, (atstart, atstop) in enumerate(ss[chain][s]):
                # see http://www.wwpdb.org/documentation/format23/sect5.html
                if s == 'H':
                    ssstr.append("HELIX  %3d %3d %3s %1s %4s%1s %3s %s %4d%s%2d%30s %5d\n"%(
                        (i+1), (i+1),
                        atstart.resn, atstart.chain, atstart.resi_number, ' ',
                        atstop.resn,  atstop.chain,  atstop.resi_number, ' ',
                        1, ' ', (atstop.resi_number -  atstart.resi_number + 1)))
                elif s == 'S':
                    ssstr.append("SHEET  %3d %3d%2d %3s %1s%4d%1s %3s %1s%4d%1s%2d\n"%(
                        (i+1), (i+1), 1, 
                        atstart.resn, atstart.chain, atstart.resi_number, '',
                        atstop.resn,  atstop.chain,  atstop.resi_number, '',   
                        0))
    return ''.join(ssstr)

def get_pdb_seqres(selection='all', quiet=1):
    '''
DESCRIPTION

    Get PDB SEQRES records for a given selection.
    '''
    prev = [None, None]
    sequences = []

    def callback(modelchain, resi, resn):
        if prev[0] != modelchain:
            prev[0] = modelchain
            sequences.append((modelchain[1], []))
        elif prev[1] == resi:
            # same residue
            return
        prev[1] = resi
        sequences[-1][1].append(resn)

    cmd.iterate('polymer & (%s)' % selection,
            'callback((model, chain), resi, resn)', space=locals())

    buf = []
    for (chain, seq) in sequences:
        numRes = len(seq)
        for i in range(0, numRes, 13):
            buf.append('SEQRES %3i %1.1s %4i  ' % ((i / 13) + 1, chain, numRes))
            for j in range(i, min(i + 13, numRes)):
                buf.append('%3.3s ' % seq[j])
            buf.append('\n')

    buf = ''.join(buf)

    if not int(quiet):
        print(buf)

    return buf

def save_pdb(filename, selection='(all)', state=-1, symm=1, ss=1, aniso=0, seqres=0, quiet=1):
    '''
DESCRIPTION

    Save the coordinates of a selection as pdb including the
    secondary structure information and, if possible, the unit
    cell. The latter requires the selction of a single object

USAGE

    save_pdb filename, selection [, state [, symm [, ss [, aniso ]]]]

ARGUMENTS

    filename = string: file path to be written

    selection = string: atoms to save {default: (all)}
                Note: to include the unit cell information you
                need to select a single object

    state = integer: state to save {default: -1 (current state)}

    symm = 0 or 1: save symmetry info if possible {default: 1}

    ss = 0 or 1: save secondary structure info {default: 1}

    aniso = 0 or 1: save ANISO records {default: 0}

SEE ALSO

    save
    '''
    _assert_package_import()
    from . import pymol_version

    selection = selector.process(selection)
    state, quiet = int(state), int(quiet)
    symm, ss     = int(symm), int(ss)
    
    filename = cmd.exp_path(filename)
    f = open(filename, 'w')
    print('REMARK 200 Generated with PyMOL and psico'.ljust(80), file=f)

    # Write sequence
    if int(seqres):
        f.write(get_pdb_seqres(selection))

    # Write the CRYST1 line if possible
    if symm:
        try:
            obj1 = cmd.get_object_list(selection)[0]
            sym = cmd.get_symmetry(obj1)
            if len(sym) != 7:
                raise
            f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%4d\n" % tuple(sym + [1]))
            if not quiet:
                print(' Info: Wrote unit cell and space group info')
        except:
            if not quiet:
                print(' Info: No crystal information')

    # Write secondary structure
    if ss:
        try:
            sss = get_pdb_sss(selection, state, quiet)
            if not sss:
                raise
            f.write(sss)
            if not quiet:
                print(' Info: Wrote secondary structure info')
        except:
            if not quiet:
                print(' Info: No secondary structure information')
    
    # Write coordinates of selection
    pdbstr = cmd.get_pdbstr(selection, state)

    # fix END records
    if state == 0 and pymol_version < 1.6:
        pdbstr = '\n'.join(line for line in pdbstr.splitlines() if line != 'END') + '\nEND\n'

    # anisotropic b-factors
    if int(aniso) and pymol_version < 1.6 and \
            cmd.get_model('first (%s)' % selection).atom[0].u_aniso[0] != 0.0:
        def mergeaniso():
            atom_it = iter(cmd.get_model(selection, state).atom)
            for line in pdbstr.splitlines(True):
                yield line
                if line[:6] in ['ATOM  ', 'HETATM']:
                    yield 'ANISOU' + line[6:28] + \
                            ''.join('%7.0f' % (u*1e4) for u in next(atom_it).u_aniso) + \
                            line[70:]
        pdbstr = ''.join(mergeaniso())

    f.write(pdbstr)
    f.close()

    if not quiet:
        print('Wrote PDB to \''+filename+'\'')

def save(filename, selection='(all)', state=-1, format='',
        ref='', ref_state=-1, quiet=1, *args, **kwargs):
    '''
ADDITIONAL NOTE

    This is an overloaded version of the 'save' command that also saves
    secondary structure and crystal information, if available.
    '''
    if not (format == 'pdb' or format == '' and filename.endswith('.pdb')):
        from pymol.exporting import save
        return save(filename, selection, state, format, ref, ref_state, quiet, *args, **kwargs)
    save_pdb(filename, selection, state, 1, 1, 0, quiet)
save.__doc__ = cmd.save.__doc__ + save.__doc__

def unittouu(string, dpi=90.0):
    '''
DESCRIPTION

    API only. Returns pixel units given a string representation of units in
    another system. Default unit is millimeter.
    '''
    import re
    uuconv = {'in':dpi, 'mm':dpi/25.4, 'cm':dpi/2.54}
    unit = 'mm'
    if isinstance(string, str) and string[-2:].isalpha():
        string, unit = string[:-2], string[-2:]
    try:
        retval = float(string)
    except:
        raise ValueError("cannot parse value from: " + str(string))
    if unit not in uuconv:
        raise ValueError("unknown unit: " + str(unit))
    return retval * uuconv[unit]

def paper_png(filename, width=100, height=0, dpi=300, ray=1):
    '''
DESCRIPTION

    Saves a PNG format image file of the current display.
    Instead of pixel dimensions, physical dimensions for
    printing (in millimeters) and DPI are specified.

USAGE

    paper_png filename [, width [, height [, dpi [, ray ]]]]

ARGUMENTS

    filename = string: filename

    width = float: width in millimeters {default: 100 = 10cm}
    width = string: width including unit (like: '10cm' or '100mm')

    height = float or string, like "width" argument. If height=0, keep
    aspect ratio {default: 0}

    dpi = float: dots-per-inch {default: 300}

    ray = 0 or 1: should ray be run first {default: 1 (yes)}

SEE ALSO

    png
    '''
    dpi, ray = float(dpi), int(ray)
    width = unittouu(width, dpi)
    height = unittouu(height, dpi)
    cmd.png(filename, width, height, dpi, ray)

def save_pdb_without_ter(filename, selection, *args, **kwargs):
    '''
DESCRIPTION

    Save PDB file without TER records. External applications like TMalign and
    DynDom stop reading PDB files at TER records, which might be undesired in
    case of missing loops.
    '''
    v = cmd.get_setting_boolean('pdb_use_ter_records')
    if v: cmd.set('pdb_use_ter_records', 0)
    cmd.save(filename, selection, *args, **kwargs)
    if v: cmd.set('pdb_use_ter_records')

## pymol command stuff

cmd.extend('save_traj', save_traj)
cmd.extend('save_pdb', save_pdb)
cmd.extend('paper_png', paper_png)

cmd.auto_arg[1].update([
    ('save_traj', cmd.auto_arg[1]['save']),
    ('save_pdb', cmd.auto_arg[1]['save']),
])

# vi:expandtab:smarttab
