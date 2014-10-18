'''
(c) 2012 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def consurfdb(code, chain='A', selection=None, palette='red_white_blue', quiet=1):
    '''
DESCRIPTION

    Color by evolutionary conservation. Writes scores to b-factor.
 
    Fetches pre-calculated conservation profile from ConSurf-DB.
    http://consurfdb.tau.ac.il/

USAGE

    consurfdb code [, chain [, selection [, palette ]]]

EXAMPLE

    fetch 1ubq, async=0
    consurfdb 3nhe, B, 1ubq

SEE ALSO

    load_consurf
    '''
    import urllib2

    code = code.upper()
    url = 'http://bental.tau.ac.il/new_ConSurfDB/DB/%s/%s/r4s.res' % (code, chain)

    try:
        handle = urllib2.urlopen(url)
    except urllib2.HTTPError:
        print ' error: no pre-calculated profile for %s/%s' % (code, chain)
        raise CmdException

    if selection is None:
        object_list = cmd.get_object_list()
        if code not in object_list:
            if code.lower() in object_list:
                code = code.lower()
            else:
                from .importing import fetch
                fetch(code, async=0)
        selection = '%s and chain %s' % (code, chain)

    load_consurf(handle, selection, palette, quiet)

def load_consurf(filename, selection, palette='red_white_blue', quiet=1):
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

    # reduced pattern that matches both r4s.res and consurf.grades
    pattern = re.compile(r'\s*(\d+)\s+([A-Y])\s+([-.0-9]+)\s')

    scores = []
    seqlist = []

    if isinstance(filename, basestring):
        handle = open(filename)
    else:
        handle = filename

    if len(cmd.get_chains(selection)) > 1:
        print ' Warning: selection spans multiple chains'

    for line in handle:
        if line.startswith('#') or line.strip() == '':
            continue
        m = pattern.match(line)
        if m is None:
            continue
        scores.append(float(m.group(3)))
        seqlist.append(m.group(2))

    selection = '(%s) and polymer' % selection
    model_ca = cmd.get_model(selection + ' and guide')
    model_seq = ''.join(one_letter.get(a.resn, 'X') for a in model_ca.atom)
    sequence = ''.join(seqlist)

    aln = needle_alignment(model_seq, sequence)
    scores_resi = dict((model_ca.atom[i].resi, scores[j])
            for (i, j) in alignment_mapping(*aln))

    cmd.alter(selection, 'b=scores.get(resi, -10)',
            space={'scores': scores_resi}, quiet=quiet)

    if palette:
        cmd.color('yellow', selection + ' and b<-9')
        if ' ' in palette:
            from .viewing import spectrumany as spectrum
        else:
            spectrum = cmd.spectrum
        spectrum('b', palette, selection + ' and b>-9.5')

cmd.extend('consurfdb', consurfdb)
cmd.extend('load_consurf', load_consurf)

# vi:expandtab:smarttab
