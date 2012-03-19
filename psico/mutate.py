'''
(c) 2010 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd

def mutate(selection, new_resn, inplace=0, sculpt=0, hydrogens='auto', mode=0, quiet=1):
    '''
DESCRIPTION

    Mutate a single residue. Does call the mutagenesis wizard non-interactively
    and tries to select the best rotamer. Can do some sculpting in the end to
    the best rotamer.

USAGE

    mutate selection, new_resn [, inplace [, sculpt [, hydrogens]]]

ARGUMENTS

    selection = string: atom selection of single residue

    new_resn = string: new residue name (3-letter or 1-letter)

    inplace = 0 or 1: work on copy of input object if 0 {default: 0}

    sculpt = 0: no sculpting {default}
    sculpt = 1: do sculpting on best rotamer
    sculpt = 2: do sculpting with neighbors

    hydrogens = string: keep, auto or none {default: auto}

    mode = 0: select rotamer with best clash score {default}
    mode = 1: take chi angles from original residue

EXAMPLE

    fetch 2x19, async=0
    select x, A/CYS`122/
    zoom x
    mutate x, LYS
    '''
    from pymol.wizard import mutagenesis
    from . import three_letter

    inplace, sculpt = int(inplace), int(sculpt)
    mode = int(mode)
    quiet = int(quiet)
    org = cmd.get_object_list(selection)[0]
    tmp = cmd.get_unused_name()
    new_resn = new_resn.upper()
    new_resn = three_letter.get(new_resn, new_resn)

    if inplace:
        cpy = org
    else:
        cpy = cmd.get_unused_name(org + '_cpy')
        cmd.create(cpy, org, -1, 1, zoom=0)

    scr = []
    cmd.iterate('first (%s)' % selection, 'scr[:] = (segi,chain,resi,resn)', space={'scr': scr})
    res = '/%s/%s/%s/%s' % tuple([cpy] + scr[:3])

    if mode == 1:
        old_resn = scr[3]
        chi_atoms = {
            'ALA': [],
            'ARG': ['C','CA','CB','CG','CD','NE','CZ'],
            'ASN': ['C','CA','CB','CG','OD1'],
            'ASP': ['C','CA','CB','CG','OD1'],
            'CYS': ['C','CA','CB','SG'],
            'GLN': ['C','CA','CB','CG','CD','OE1'],
            'GLU': ['C','CA','CB','CG','CD','OE1'],
            'GLY': [],
            'HIS': ['C','CA','CB','CG','ND1'],
            'ILE': ['C','CA','CB','CG1','CD1'],
            'LEU': ['C','CA','CB','CG','CD1'],
            'LYS': ['C','CA','CB','CG','CD','CE','NZ'],
            'MET': ['C','CA','CB','CG','SD','CE'],
            'MSE': ['C','CA','CB','CG','SE','CE'],
            'PHE': ['C','CA','CB','CG','CD1'],
            'PRO': [],
            'SER': ['C','CA','CB','OG'],
            'THR': ['C','CA','CB','OG1'],
            'TRP': ['C','CA','CB','CG','CD2'],
            'TYR': ['C','CA','CB','CG','CD1'],
            'VAL': ['C','CA','CB','CG2'],
        }
        atoms = [res + '/' + name for name in chi_atoms.get(old_resn, [])]
        old_chi = [cmd.get_dihedral(*args) for args in zip(atoms, atoms[1:], atoms[2:], atoms[3:])]

    cmd.remove('%s and not name CA+C+N+O+OXT' % (res))

    # start the wizard to count the number of rotamers for this residue
    cmd.wizard("mutagenesis")
    cmd.get_wizard().set_mode(new_resn)
    cmd.get_wizard().set_hyd(hydrogens)
    cmd.get_wizard().do_select("("+res+")")

    def get_best_state_bump():
        best_state = (1, 1e9)
        cmd.create(tmp, '%s and not name CA+C+N+O or (%s within 8.0 of (%s and name CB))' % \
                (mutagenesis.obj_name, cpy, mutagenesis.obj_name), zoom=0, singletons=1)
        cmd.bond('name CB and %s in %s' % (tmp, mutagenesis.obj_name),
                'name CA and %s in %s' % (tmp, res))
        cmd.sculpt_activate(tmp)
        for i in range(1, cmd.count_states(tmp)+1):
            score = cmd.sculpt_iterate(tmp, state=i)
            if not quiet:
                print 'Frame %d Score %.2f' % (i, score)
            if score < best_state[1]:
                best_state = (i, score)
        cmd.delete(tmp)
        if not quiet:
            print ' Best: Frame %d Score %.2f' % best_state
        return best_state

    if cmd.count_states(mutagenesis.obj_name) < 2 or mode == 1:
        best_state = (1, -1.0)
    else:
        best_state = get_best_state_bump()
    cmd.frame(best_state[0])
    cmd.get_wizard().apply()
    cmd.wizard()

    if mode == 1:
        atoms = [res + '/' + name for name in chi_atoms.get(new_resn, [])]
        for args in zip(atoms, atoms[1:], atoms[2:], atoms[3:], old_chi):
            cmd.set_dihedral(*args)
        cmd.unpick()

    if sculpt > 0:
        residue_sculpting(res, cpy)
    
    return cpy

def mutate_all(selection, new_resn, inplace=1, sculpt=0, *args, **kwargs):
    '''
DESCRIPTION

    Mutate all residues in selection. By default do mutation in-place (unlike
    the 'mutate' command which by default works on a copy).

    FOR SCULPTING ONLY SUPPORTS SELECTIONS WITHIN THE SAME MODEL!

SEE ALSO

    mutate
    '''
    inplace, sculpt = int(inplace), int(sculpt)

    if sculpt and len(cmd.get_object_list('(' + selection + ')')) > 1:
        print ' Error: Sculpting in multiple models not supported'
        raise CmdException

    kwargs.pop('_self', None)
    sele_list = set()
    cmd.iterate(selection,
            'sele_list.add("/%s/%s/%s/%s" % (model, segi, chain, resi))',
            space={'sele_list': sele_list})
    for sele in sele_list:
        mutate(sele, new_resn, inplace, sculpt and not inplace, *args, **kwargs)
    if sculpt and inplace:
        residue_sculpting('(' + ' '.join(sele_list) + ')')

def residue_sculpting(sele, model=''):
    '''
    Relax the given selection.

    SO FAR ONLY SUPPORTS SELECTIONS WITHIN THE SAME MODEL!

    Do 100 iterations, 75 of them with all terms but low VDW weights,
    and 25 with only local geometry terms. With default VDW weights and
    atom clashes, the structure gets distorted very easily!
    '''
    from pymol import selector
    sele = selector.process(sele)
    org = cmd.get_object_list(sele)[0]
    if len(model) == 0:
        model = org
    elif model != org:
        print 'org -> model:', org, '->', model
        sele = sele.replace('(%s)' % org, '(%s)' % model)
    cmd.protect('(not %s) or name CA+C+N+O+OXT' % (sele))
    cmd.sculpt_activate(model)
    cmd.set('sculpt_vdw_weight', 0.25, model) # Low VDW forces
    cmd.set('sculpt_field_mask', 0x1FF, model) # Default
    if sculpt == 2:
        cmd.sculpt_iterate(model, cycles=25)
        cmd.deprotect('(byres (%s within 6.0 of %s)) and (not name CA+C+N+O+OXT)' % \
                (model, sele))
        cmd.sculpt_iterate(model, cycles=50)
    else:
        cmd.sculpt_iterate(model, cycles=75)
    cmd.set('sculpt_field_mask', 0x01F, model) # Local Geometry Only
    cmd.sculpt_iterate(model, cycles=25)
    cmd.unset('sculpt_vdw_weight', model)
    cmd.unset('sculpt_field_mask', model)
    cmd.sculpt_deactivate(model)
    cmd.deprotect()

cmd.extend('mutate', mutate)
cmd.extend('mutate_all', mutate_all)

cmd.auto_arg[0].update({
    'mutate'         : cmd.auto_arg[0]['align'],
    'mutate_all'     : cmd.auto_arg[0]['align'],
})

# vi: expandtab:smarttab
