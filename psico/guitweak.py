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

def colorramps(self_cmd, sele):
    '''Menu for coloring atoms and surfaces by color ramps'''
    return [[ 1, name, 'cmd.color("%s", "%s");'
        'cmd.set("surface_color", "%s", "%s")' % (name, sele, name, sele) ]
        for name in self_cmd.get_names_of_type('object:')]

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
def ramp_action(self_cmd, sele):
    r = [
        [ 0, '', '' ],
        [ 1, 'levels', [
            [ 1, 'Range +/- %.1f' % (L),
                'psico.creating.ramp_levels("%s", [%f, 0, %f])' % (sele, -L, L) ]
            for L in [1, 2, 5, 10, 20, 50, 100]
        ]],
    ]
    return r

@menuappend
def all_colors(self_cmd, sele):
    ramps = colorramps(self_cmd, sele)
    if len(ramps) == 0:
        return []
    return [[ 0, '', '' ], [ 1, 'ramps', ramps ]]

@menuappend
def slice_color(self_cmd, sele):
    return colorramps(self_cmd, sele)

@menuappend
def map_volume(self_cmd, sele):
    m = 'psico.electrostatics'
    return [
        [ 1, 'electrostatics', m + '.volume_esp("'+sele+'_volume_esp","'+sele+'")' ],
        [ 1, 'difference density', m + '.volume_fofc("'+sele+'_volume_fofc","'+sele+'")' ],
    ]

# vi: ts=4:sw=4:smarttab:expandtab
