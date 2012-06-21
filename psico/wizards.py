'''
Interactive wizards

(c) 2011-2012 Thomas Holder

License: BSD-2-Clause
'''

import sys
from pymol import cmd
from pymol.wizard import Wizard

class Sspick(Wizard):
    '''
    Secondary structure element picker
    '''

    def __init__(self, _self=cmd):

        cmd.unpick();
        Wizard.__init__(self, _self)

        self.mouse_selection_mode = cmd.get_setting_legacy('mouse_selection_mode')
        cmd.set('mouse_selection_mode',0) # set selection mode to atomic
        cmd.deselect() # disable the active selection (if any)

        self.error = None
        self.selection_mode = 1 if self.mouse_selection_mode == 6 else 0
        self.selection_modes = [
            'Residues',
            'C-alphas',
        ]

        smm = []
        smm.append([ 2, 'Selection Mode', '' ])
        for i, label in enumerate(self.selection_modes):
            smm.append([1, label, 'cmd.get_wizard().set_mode(%d)' % i])
        self.menu['selection_mode'] = smm

        self.name = None

    def set_mode(self, i):
        self.selection_mode = i
        cmd.refresh_wizard()

    def get_panel(self):
        return [
            [1, 'SS Picker', ''],
            [3, self.selection_modes[self.selection_mode], 'selection_mode'],
            [2, 'Done', 'cmd.set_wizard()'],
            ]

    def cleanup(self):
        cmd.set('mouse_selection_mode', self.mouse_selection_mode)

    def get_prompt(self):
        self.prompt = ['Pick an atom']
        if self.error is not None:
            self.prompt.append(self.error)
        return self.prompt

    def do_select(self, name): # map selects into picks
        from .selecting import select_sspick
        if self.name not in cmd.get_names('selections', enabled_only=1):
            self.name = cmd.get_unused_name('ss')
        select_sspick(name, self.name, self.selection_mode)
        cmd.enable(self.name)
        cmd.refresh_wizard()

    def do_pick(self, bondFlag):
        self.do_select('(pk1)')
        cmd.unpick()

# enabling wizards with the "wizard" command is a bit ugly...

names = ['sspick']

for name in names:
    sys.modules['pymol.wizard.' + name] = sys.modules[__name__]

# tab-completion for "wizard"

def names_sc():
    import os, pymol.wizard
    names_glob = [name[:-3] for p in pymol.wizard.__path__
            for name in os.listdir(p) if name.endswith('.py')]
    return cmd.Shortcut(names_glob + names)

cmd.auto_arg[0]['wizard'] = [ names_sc, 'wizard', '' ]

# vi:expandtab:smarttab:sw=4
