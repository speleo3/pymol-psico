'''
(c) 2010-2012 Thomas Holder

License: BSD-2-Clause

Extensions to the PyMOL GUI like menu items
'''

from pymol import menu

# originial (not overloaded) menu functions
orig = dict((name, getattr(menu, name)) for name in [
    'mol_generate',
    'ramp_action',
    ])

def mol_generate(self_cmd, sele):
    try:
        from epymol import rigimol
        cmd = 'morpheasy'
    except ImportError:
        cmd = 'morpheasy_linear'
    r = orig['mol_generate'](self_cmd, sele) + [
        [ 0, '', '' ],
        [ 1, 'electrostatics (APBS)', 'psico.electrostatics.apbs_surface("'+sele+'")' ],
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
        [ 0, '', '' ],
        [ 1, 'morphing to ...', [
            [ 1, a, 'psico.morphing.%s(%s,%s)' % (cmd, repr(sele), repr(a)) ]
                for a in self_cmd.get_object_list()[0:25] if a != sele
        ]],
    ]
    return r

def ramp_action(self_cmd, sele):
    r = orig['ramp_action'](self_cmd, sele) + [
        [ 0, '', '' ],
        [ 1, 'levels', [
            [ 1, 'Range +/- %.1f' % (L),
                'psico.creating.ramp_levels("%s", [%f, 0, %f])' % (sele, -L, L) ]
            for L in [1, 2, 5, 10, 20, 50, 100]
        ]],
    ]
    return r

# overload menu functions
for name in orig:
    setattr(menu, name, globals()[name])

# vi: ts=4:sw=4:smarttab:expandtab
