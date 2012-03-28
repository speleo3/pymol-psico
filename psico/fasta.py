'''
(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd

def fasta(selection='(all)', gapped=1, wrap=70):
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
    selection = '(%s) and polymer' % (selection)
    for obj in cmd.get_object_list(selection):
        for chain in cmd.get_chains('%s and (%s)' % (obj, selection)):
            seq = []
            model = cmd.get_model('/%s//%s//CA and (%s)' % (obj, chain, selection))
            prev_resi = 999999999
            for atom in model.atom:
                if gapped:
                    gap_len = max(0, atom.resi_number - prev_resi - 1)
                    seq.extend('-' * gap_len)
                    prev_resi = atom.resi_number
                seq.append(one_letter.get(atom.resn, 'X'))
            print '>%s_%s' % (obj, chain)
            if wrap < 1:
                print ''.join(seq)
                continue
            for i in range(0, len(seq), wrap):
                print ''.join(seq[i:i+wrap])

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
        model = cmd.get_model('/%s////CA and polymer and (%s)' % (obj, selection))
        for atom in model.atom:
            if prev_coord is not None and cpv.distance(atom.coord, prev_coord) > 4.0:
                seq.append('/\n')
            prev_coord = atom.coord
            seq.append(one_letter.get(atom.resn, 'X'))
        seq.append('*')
        print '>P1;%s' % (obj)
        print 'structure:%s:%s:%s:%s:%s::::' % (obj,
                model.atom[0].resi,model.atom[0].chain,
                model.atom[-1].resi,model.atom[-1].chain)
        if wrap < 1:
            print ''.join(seq)
            continue
        for i in range(0, len(seq), wrap):
            print ''.join(seq[i:i+wrap])

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
        color = '#%02x%02x%02x' % (c[0]*255, c[1]*255, c[2]*255)
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
    print >> handle, '<html><body style="font-family:monospace">'
    print >> handle, ''.join(html)
    print >> handle, '</body></html>'
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
