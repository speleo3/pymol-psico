'''
(c) 2010-2012 Thomas Holder

License: BSD-2-Clause

Extensions to the PyMOL GUI like menu items
'''

from pymol import menu

dummy = lambda c, s: []

def menuappend(f):
    '''Decorator for overloading menu functions by appending to them'''
    orig = getattr(menu, f.__name__, dummy)
    wrapped = lambda c, s: orig(c, s) + f(c, s)
    setattr(menu, f.__name__, wrapped)
    return wrapped

@menuappend
def mol_generate(self_cmd, sele):
    r = [
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
    ]
    return r

@menuappend
def map_volume(self_cmd, sele):
    m = 'psico.electrostatics'
    return [
        [ 0, '', '' ],
        [ 1, 'electrostatics', m + '.volume_esp("'+sele+'_volume_esp","'+sele+'")' ],
        [ 1, 'difference density', m + '.volume_fofc("'+sele+'_volume_fofc","'+sele+'")' ],
    ]

@menuappend
def presets(self_cmd, sele):
    return [
        [ 0, '', '' ],
        [ 1, 'hydropathy', 'psico.aaindex.hydropathy2b("polymer & (' + sele + ')",quiet=0)' ],
    ]

# vi: ts=4:sw=4:smarttab:expandtab
