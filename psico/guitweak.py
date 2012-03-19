'''
(c) 2010 Thomas Holder

License: BSD-2-Clause

Extensions to the PyMOL GUI like menu items
'''

from pymol import menu

x__mol_generate = menu.mol_generate

def mol_generate(self_cmd, sele):
    r = x__mol_generate(self_cmd, sele) + [
        [ 0, '', '' ],
        [ 1, 'biological unit', 'psico.xtal.biomolecule("'+sele+'")' ],
        [ 1, 'supercell with mates', [
            [ 2, 'Supercell:', '' ],
            [ 1, '1x1x1', 'psico.xtal.supercell(1,1,1,"'+sele+'")' ],
            [ 1, '2x2x2', 'psico.xtal.supercell(2,2,2,"'+sele+'")' ],
            [ 1, '2x1x1', 'psico.xtal.supercell(2,1,1,"'+sele+'")' ],
            [ 1, '1x2x1', 'psico.xtal.supercell(1,2,1,"'+sele+'")' ],
            [ 1, '1x1x2', 'psico.xtal.supercell(1,1,2,"'+sele+'")' ],
        ]],
        [ 1, 'sequence', [
            [ 2, 'Sequence:', '' ],
            [ 1, 'Fasta', 'psico.fasta.fasta("'+sele+'")' ],
            [ 1, 'PIR', 'psico.fasta.pir("'+sele+'")' ],
        ]],
    ]
    try:
        from epymol import rigimol
        morph_submenu = []
        for a in self_cmd.get_object_list()[0:25]:
            if a != sele:
                morph_submenu.append([ 1, a, 'psico.morpheasy.morpheasy(%s,%s)' % (repr(sele), repr(a)) ])
        r.extend([
            [ 0, '', '' ],
            [ 1, 'morphing to ...', morph_submenu ]
        ])
    except ImportError:
        pass
    return r

menu.mol_generate = mol_generate

# vi: ts=4:sw=4:smarttab:expandtab
