'''
(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

if not __name__.endswith('.fasta'):
    raise Exception("Must do 'import psico.fasta' instead of 'run ...'")

import sys

from pymol import cmd

def fasta(selection='(all)', gapped=1, wrap=70, filename='', quiet=1):
    '''
DESCRIPTION

    Print sequence in FASTA format

ARGUMENTS

    selection = string: atom selection {default: all}

    gapped = integer: put missing residues as dashes into the sequence {default: 1}

    wrap = integer: wrap lines, 0 for no wrapping {default: 70}

SEE ALSO

    pir, pymol.exporting.get_fastastr
    '''
    from . import one_letter
    gapped, wrap = int(gapped), int(wrap)

    class prev:
        key = None
        col = 0

    if filename:
        out = open(filename, 'w')
    else:
        out = sys.stdout

    def write(c):
        if wrap > 0 and prev.col == wrap:
            out.write('\n')
            prev.col = 0
        prev.col += 1
        out.write(c)

    def callback(key, resv, resn):
        if key != prev.key:
            # different chain
            if prev.col > 0:
                out.write('\n')
            out.write('>%s_%s\n' % key)
            prev.key = key
            prev.col = 0
        elif resv == prev.resv:
            # same residue
            return
        elif gapped:
            for _ in range(resv - prev.resv - 1):
                # gap
                write('-')
        prev.resv = resv
        write(one_letter.get(resn, 'X'))

    cmd.iterate('(%s) & polymer' % (selection),
            '_cb((model, chain), resv, resn)',
            space={'_cb': callback})

    if prev.col > 0:
        out.write('\n')

    if filename:
        if not int(quiet):
            print(' Wrote sequence to "%s"' % filename)
        out.close()

def pir(selection='(all)', wrap=70):
    '''
DESCRIPTION

    Print sequence in PIR format

SEE ALSO

    fasta
    '''
    from . import one_letter
    from chempy import cpv
    wrap = int(wrap)
    for obj in cmd.get_object_list('(' + selection + ')'):
        seq = []
        prev_coord = None
        model = cmd.get_model('/%s////CA and guide and (%s)' % (obj, selection))
        for atom in model.atom:
            if prev_coord is not None and cpv.distance(atom.coord, prev_coord) > 4.0:
                seq.append('/\n')
            prev_coord = atom.coord
            seq.append(one_letter.get(atom.resn, 'X'))
        seq.append('*')
        print('>P1;%s' % (obj))
        print('structure:%s:%s:%s:%s:%s::::' % (obj,
                model.atom[0].resi,model.atom[0].chain,
                model.atom[-1].resi,model.atom[-1].chain))
        if wrap < 1:
            print(''.join(seq))
            continue
        for i in range(0, len(seq), wrap):
            print(''.join(seq[i:i+wrap]))

def save_colored_fasta(filename, selection='(all)', gapped=1, quiet=1):
    '''
DESCRIPTION

    Save a html file with colored (by C-alpha atoms) fasta sequence.
    '''
    from . import one_letter
    from pymol import Scratch_Storage
    gapped = int(gapped)
    selection = '(%s) and polymer and guide' % (selection)
    html = []
    stored = Scratch_Storage()
    def callback(resv, resn, color):
        if stored.resv is None:
            stored.resv = resv - (resv % 70)
        if gapped:
            while stored.resv+1 < resv:
                callback(stored.resv+1, '-', 25)
        stored.resv += 1
        if stored.resv % 70 == 1:
            html.append(('</font>\n<br>%4d <font>' % (resv)).replace(' ', '&nbsp;'))
            stored.color = None
        c = cmd.get_color_tuple(color)
        color = '#%02x%02x%02x' % tuple(int(0xFF * v) for v in c)
        aa = one_letter.get(resn, '-')
        if color != stored.color:
            html.append('</font><font color="' + color + '">')
            stored.color = color
        html.append(aa)
    for obj in cmd.get_object_list('(' + selection + ')'):
        for chain in cmd.get_chains('model %s and (%s)' % (obj, selection)):
            sele = 'model %s and chain "%s" and (%s)' % (obj, chain, selection)
            html.append('\n<br>&gt;%s_%s<font>' % (obj, chain))
            stored.resv = None if gapped else 0
            stored.color = None
            cmd.iterate(sele, 'callback(resv, resn, color)', space=locals())
            html.append('</font>')
    handle = open(filename, 'w')
    print('<html><body style="font-family:monospace">', file=handle)
    print(''.join(html), file=handle)
    print('</body></html>', file=handle)
    handle.close()

cmd.extend('fasta', fasta)
cmd.extend('pir', pir)
cmd.extend('save_colored_fasta', save_colored_fasta)

cmd.auto_arg[0].update({
    'fasta'          : cmd.auto_arg[0]['zoom'],
    'pir'            : cmd.auto_arg[0]['zoom'],
})
cmd.auto_arg[1].update([
    ('save_colored_fasta', cmd.auto_arg[0]['zoom']),
])

# vi: ts=4:sw=4:smarttab:expandtab
