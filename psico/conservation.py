'''
(c) 2012 Thomas Holder

License: BSD-2-Clause
'''

if not __name__.endswith('.conservation'):
    raise Exception("Must do 'import psico.conservation' instead of 'run ...'")

from pymol import cmd, CmdException


def consurfdb(code, chain='A', selection=None, palette='red_white_blue', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Color by evolutionary conservation. Writes scores to b-factor.

    Fetches pre-calculated conservation profile from ConSurf-DB.
    http://consurfdb.tau.ac.il/

USAGE

    consurfdb code [, chain [, selection [, palette ]]]

EXAMPLE

    fetch 1ubq, async=0
    consurfdb 5nvg, A, 1ubq

SEE ALSO

    load_consurf
    '''
    import ssl
    import urllib.request as urllib2

    code = code.upper()
    url = f'https://consurfdb.tau.ac.il/DB/{code}{chain}/{code}{chain}_consurf_summary.txt'
    context = ssl._create_unverified_context()

    try:
        handle = urllib2.urlopen(url, context=context)
    except urllib2.HTTPError as ex:
        raise CmdException(
            f'no pre-calculated profile for {code}/{chain}') from ex

    if selection is None:
        object_list = _self.get_object_list()
        if code not in object_list:
            if code.lower() in object_list:
                code = code.lower()
            else:
                from .importing import fetch
                fetch(code, _self=_self)
        selection = '%s and chain %s' % (code, chain)

    load_consurf(handle, selection, palette, quiet, _self=_self)


def load_consurf(filename, selection, palette='red_white_blue', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Color by evolutionary conservation. Writes scores to b-factor.

    You need a "r4s.res" or "consurf.grades" input file.

USAGE

    load_consurf filename, selection [, palette ]

SEE ALSO

    consurfdb
    '''
    import re
    from .seqalign import needle_alignment, alignment_mapping
    from . import one_letter

    # reduced pattern that matches r4s.res, consurf.grades and consurf_summary.txt
    pattern = re.compile(r'\s*(\d+)\s+([A-Y])\s+(?:-\s+|\S+:\w*\s+)?([-.0-9]+)\s')

    scores = []
    seqlist = []

    if isinstance(filename, (str, bytes)):
        handle = open(filename)
    else:
        handle = filename

    if len(_self.get_chains(selection)) > 1:
        print(' Warning: selection spans multiple chains')

    for line in handle:
        if not isinstance(line, str):
            line = line.decode()
        if line.startswith('#') or line.strip() == '':
            continue
        m = pattern.match(line)
        if m is None:
            continue
        scores.append(float(m.group(3)))
        seqlist.append(m.group(2))

    selection = '(%s) and polymer' % selection
    model_ca = _self.get_model(selection + ' and guide')
    model_seq = ''.join(one_letter.get(a.resn, 'X') for a in model_ca.atom)
    sequence = ''.join(seqlist)

    aln = needle_alignment(model_seq, sequence)
    scores_resi = dict((model_ca.atom[i].resi, scores[j])
            for (i, j) in alignment_mapping(*aln))

    _self.alter(selection, 'b=scores.get(resi, -10)',
            space={'scores': scores_resi}, quiet=quiet)

    if palette:
        _self.color('yellow', selection + ' and b<-9')
        _self.spectrum('b', palette, selection + ' and b>-9.5')


cmd.extend('consurfdb', consurfdb)
cmd.extend('load_consurf', load_consurf)

# vi:expandtab:smarttab
