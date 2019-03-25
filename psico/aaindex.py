'''
Python parser for AAindex: Amino Acid Index Database

http://www.genome.jp/aaindex/

(c) 2010-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from __future__ import print_function

from pymol import cmd, CmdException

def _assert_package_import():
    if not __name__.endswith('.aaindex'):
        raise CmdException("Must do 'import psico.aaindex' instead of 'run ...'")

_aaindex1 = None
def _get_aaindex1():
    '''
    Get a global AAindex instance for "fetch_path/aaindex1"
    '''
    global _aaindex1
    if _aaindex1 is None:
        _aaindex1 = AAindex(1, cmd.get('fetch_path'), 0)
    return _aaindex1

class AAindex(dict):
    '''
    Represents an aaindex database file as key => record mapping.
    '''

    def __init__(self, index='', path=None, quiet=1):
        '''
        Read in the aaindex files. You need to run this (once) before you can
        access any records. If the files are not within the current directory,
        you need to specify the correct directory path. By default aaindex1 is
        read in.
        '''
        import os

        if not index:
            return

        index = str(index)

        if path is None:
            for path in [os.path.dirname(__file__), '.']:
                if os.path.exists(os.path.join(path, 'aaindex1')):
                    break

        parse = lambda a,b: self.parse(os.path.join(path, a), b, quiet)

        if '1' in index:
            parse('aaindex1', Record)
        for i in '23':
            if i in index:
                parse('aaindex' + i, MatrixRecord)

    def search(self, pattern, searchtitle=False, casesensitive=False):
        '''
        Search for pattern in description and title (optional) of all records and
        return matched records as list. By default search case insensitive.
        '''
        if casesensitive:
            whatcase = lambda i: i
        else:
            pattern = pattern.lower()
            whatcase = lambda i: i.lower()
        matches = []
        for record in self.values():
            if pattern in whatcase(record.desc) or searchtitle and pattern in whatcase(record.title):
                matches.append(record)
        return matches

    def grep(self, pattern):
        '''
        Search for pattern in title and description of all records (case
        insensitive) and print results on standard output.
        '''
        for record in self.search(pattern):
            print(record)

    def parse(self, filename, rec, quiet=True):
        '''
        Parse aaindex input file. "rec" must be "Record" for aaindex1 and
        "MarixRecord" for aaindex2 and aaindex3.
        '''
        import os

        def _float_or_None(x):
            if x == 'NA' or x == '-':
                return None
            return float(x)

        if not os.path.exists(filename):
            try:
                import urllib.request as urllib
            except ImportError:
                import urllib
            url = 'ftp://ftp.genome.jp/pub/db/community/aaindex/' + os.path.basename(filename)
            if not quiet:
                print('Downloading "%s"' % (url))
            filename = urllib.urlretrieve(url, filename)[0]
            if not quiet:
                print('Saved to "%s"' % (filename))

        current = rec()
        lastkey = None

        for line in open(filename):
            key = line[0:2]
            if key[0] == ' ':
                key = lastkey

            if key == '//':
                self[current.key] = current
                current = rec()
            elif key == 'H ':
                current.key = line[2:].strip()
            elif key == 'R ':
                current.ref += line[2:]
            elif key == 'D ':
                current.desc += line[2:]
            elif key == 'A ':
                current.authors += line[2:]
            elif key == 'T ':
                current.title += line[2:]
            elif key == 'J ':
                current.journal += line[2:]
            elif key == '* ':
                current.comment += line[2:]
            elif key == 'C ':
                a = line[2:].split()
                for i in range(0, len(a), 2):
                    current.correlated[a[i]] = float(a[i+1])
            elif key == 'I ':
                a = line[1:].split()
                if a[0] != 'A/L':
                    current.extend(map(_float_or_None, a))
                elif list(Record.aakeys) != [i[0] for i in a] + [i[-1] for i in a]:
                    print('Warning: wrong amino acid sequence for', current.key)
                else:
                    try:
                        assert list(Record.aakeys[:10]) == [i[0] for i in a]
                        assert list(Record.aakeys[10:]) == [i[2] for i in a]
                    except:
                        print('Warning: wrong amino acid sequence for', current.key)
            elif key =='M ':
                a = line[2:].split()
                if a[0] == 'rows':
                    if a[4] == 'rows':
                        a.pop(4)
                    assert a[3] == 'cols' and len(a) == 6
                    i = 0
                    for aa in a[2]:
                        current.rows[aa] = i
                        i += 1
                    i = 0
                    for aa in a[5]:
                        current.cols[aa] = i
                        i += 1
                else:
                    current.extend(map(_float_or_None, a))
            elif len(key.rstrip('\r\n')) < 2:
                pass
            elif not quiet:
                print('Warning: line starts with "%s"' % (key))

            lastkey = key
    
    def __repr__(self):
        return '<%s %s>' % (self.__class__.__name__, ' '.join(self))

class Record(object):
    '''
    Amino acid index (AAindex) Record
    '''
    aakeys = 'ARNDCQEGHILKMFPSTWYV'
    def __init__(self):
        self.key = None
        self.desc = ''
        self.ref = ''
        self.authors = ''
        self.title = ''
        self.journal = ''
        self.correlated = dict()
        self.index = dict()
        self.comment = ''
    def extend(self, row):
        i = len(self.index)
        for x in row:
            self.index[self.aakeys[i]] = x
            i += 1
    def get(self, aai, aaj=None, d=None):
        assert aaj is None
        return self.index.get(aai, d)
    def __getitem__(self, aai):
        return self.get(aai)
    def median(self):
        x = sorted([_f for _f in list(self.index.values()) if _f])
        half = len(x) // 2
        if len(x) % 2 == 1:
            return x[half]
        return (x[half-1] + x[half])/2.0
    def __repr__(self):
        return '<%s %s>' % (self.__class__.__name__, self.key)
    def __str__(self):
        desc = self.desc.replace('\n', ' ').strip()
        return '%s(%s: %s)' % (self.__class__.__name__, self.key, desc)

class MatrixRecord(Record):
    '''
    Matrix record for mutation matrices or pair-wise contact potentials
    '''
    def __init__(self):
        Record.__init__(self)
        self.index = []
        self.rows = dict()
        self.cols = dict()
    def extend(self, row):
        self.index.append(row)
    def _get(self, aai, aaj):
        i = self.rows[aai]
        j = self.cols[aaj]
        return self.index[i][j]
    def get(self, aai, aaj, d=None):
        try:
            return self._get(aai, aaj)
        except:
            pass
        try:
            return self._get(aaj, aai)
        except:
            return d
    def __getitem__(self, aaij):
        return self.get(aaij[0], aaij[1])
    def median(self):
        x = []
        for y in self.index:
            x.extend([_f for _f in y if _f])
        x.sort()
        half = len(x) // 2
        if len(x) % 2 == 1:
            return x[half]
        return (x[half - 1] + x[half]) / 2.0

def aaindex2b(key, selection='all', var='b', quiet=1):
    '''
DESCRIPTION

    Looks up the Amino Acid Index from http://www.genome.jp/aaindex/
    for the given key and assignes b-factors to the given selection. Unknown
    residues get the average index value assigned.

USAGE

    aaindex2b key [, selection ]

ARGUMENTS

    key = string: Key of AAindex entry

    selection = string: atoms to assign b-factors {default: (all)}

EXAMPLE

    # Hydropathy index by Kyte-Doolittle
    aaindex2b KYTJ820101
    spectrumany b, white yellow forest
    show surface

SEE ALSO

    hydropathy2b
    '''
    _assert_package_import()
    from . import one_letter

    quiet = int(quiet)

    aaindex = _get_aaindex1()

    try:
        entry = aaindex[key]
    except KeyError:
        print(' Error: No such key in AAindex:', key)
        raise CmdException

    median = entry.median()

    if not quiet:
        print(entry.desc.strip())

    def lookup(resn):
        aa = one_letter.get(resn, 'X')
        value = entry.get(aa)
        if value is None:
            return median
        return value

    cmd.alter(selection, var + '=lookup(resn)', space=locals())

def hydropathy2b(selection='all', var='b', quiet=1):
    '''
DESCRIPTION

    Load Kyte-Doolittle hydropathy index values into b-factor column.

SEE ALSO

    aaindex2b
    '''
    r = aaindex2b('KYTJ820101', selection, var, quiet)

    if not quiet:
        from .viewing import spectrumany
        spectrumany(var, 'white forest', selection)

    return r

# commands and tab-completion of arguments

def _aaindexkey_sc():
    aaindex = _get_aaindex1()
    return cmd.Shortcut(aaindex.keys())

cmd.extend('aaindex2b', aaindex2b)
cmd.extend('hydropathy2b', hydropathy2b)

cmd.auto_arg[0].update([
    ('aaindex2b', [ _aaindexkey_sc, 'aaindexkey', ', ' ]),
    ('hydropathy2b', cmd.auto_arg[1]['select']),
])
cmd.auto_arg[1].update([
    ('aaindex2b', cmd.auto_arg[1]['select']),
])

# vi: ts=4:sw=4:smarttab:expandtab
