'''
(c) 2010-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

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
    mode = 2: first rotamer

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
        old_chi = []
        for args in zip(atoms, atoms[1:], atoms[2:], atoms[3:]):
            try:
                old_chi.append(cmd.get_dihedral(*args))
            except:
                break

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
                print('Frame %d Score %.2f' % (i, score))
            if score < best_state[1]:
                best_state = (i, score)
        cmd.delete(tmp)
        if not quiet:
            print(' Best: Frame %d Score %.2f' % best_state)
        return best_state

    if cmd.count_states(mutagenesis.obj_name) < 2 or mode > 0:
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
        sculpt_relax(res, 0, sculpt == 2, cpy)

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
        print(' Error: Sculpting in multiple models not supported')
        raise CmdException

    kwargs.pop('_self', None)
    sele_list = set()
    cmd.iterate(selection,
            'sele_list.add("/%s/%s/%s/%s" % (model, segi, chain, resi))',
            space={'sele_list': sele_list})
    for sele in sele_list:
        mutate(sele, new_resn, inplace, sculpt and not inplace, *args, **kwargs)
    if sculpt and inplace:
        sculpt_relax('(' + ' '.join(sele_list) + ')', 0, sculpt == 2)

def sculpt_relax(selection, backbone=1, neighbors=0, model=None, cycles=100,
        state=0, quiet=1):
    '''
DESCRIPTION

    Relax the given selection.

    SO FAR ONLY SUPPORTS SELECTIONS WITHIN THE SAME MODEL!

    Do 100 iterations, 75 of them with all terms but low VDW weights,
    and 25 with only local geometry terms. With default VDW weights and
    atom clashes, the structure gets distorted very easily!

USAGE

    sculpt_relax selection [, backbone [, neighbors [, model [, cycles ]]]]
    '''
    from pymol import selector

    backbone, neighbors = int(backbone), int(neighbors)
    cycles, state, quiet = int(cycles), int(state), int(quiet)

    sele = selector.process(selection)
    org = cmd.get_object_list(sele)[0]
    if model is None:
        model = org
    elif model != org:
        sele = sele.replace('(%s)' % org, '(%s)' % model)

    cmd.protect()
    cmd.deprotect(sele)
    if not backbone:
        cmd.protect('name CA+C+N+O+OXT')

    cmd.sculpt_activate(model, state)
    cmd.set('sculpt_vdw_weight', 0.25, model) # Low VDW forces
    cmd.set('sculpt_field_mask', 0x1FF, model) # Default

    if neighbors:
        cmd.sculpt_iterate(model, state, int(cycles * 0.25))
        cmd.deprotect('byres (%s within 6.0 of (%s))' % (model, sele))
        if not backbone:
            cmd.protect('name CA+C+N+O+OXT')
        cmd.sculpt_iterate(model, state, cycles=int(cycles * 0.50))
    else:
        cmd.sculpt_iterate(model, state, int(cycles * 0.75))

    cmd.set('sculpt_field_mask', 0x01F, model) # Local Geometry Only
    cmd.sculpt_iterate(model, state, int(cycles * 0.25))

    cmd.unset('sculpt_vdw_weight', model)
    cmd.unset('sculpt_field_mask', model)
    cmd.sculpt_deactivate(model)
    cmd.deprotect()

def add_missing_atoms(selection='all', cycles=200, quiet=1):
    '''
DESCRIPTION

    Mutate those residues to themselves which have missing atoms

SEE ALSO

    stub2ala
    '''
    from collections import defaultdict
    from chempy import fragments

    cycles, quiet = int(cycles), int(quiet)

    reference = {
        'ALA': set(['CB']),
        'ARG': set(['CB', 'CG', 'NE', 'CZ', 'NH1', 'NH2', 'CD']),
        'ASN': set(['CB', 'CG', 'OD1', 'ND2']),
        'ASP': set(['CB', 'CG', 'OD1', 'OD2']),
        'CYS': set(['CB', 'SG']),
        'GLN': set(['CB', 'CG', 'CD', 'NE2', 'OE1']),
        'GLU': set(['CB', 'CG', 'OE2', 'CD', 'OE1']),
        'GLY': set([]),
        'HIS': set(['CE1', 'CB', 'CG', 'CD2', 'ND1', 'NE2']),
        'ILE': set(['CB', 'CD1', 'CG1', 'CG2']),
        'LEU': set(['CB', 'CG', 'CD1', 'CD2']),
        'LYS': set(['CB', 'CG', 'NZ', 'CE', 'CD']),
        'MET': set(['CB', 'CG', 'CE', 'SD']),
        'PHE': set(['CE1', 'CB', 'CG', 'CZ', 'CD1', 'CD2', 'CE2']),
        'PRO': set(['CB', 'CG', 'CD']),
        'SER': set(['OG', 'CB']),
        'THR': set(['CB', 'OG1', 'CG2']),
        'TRP': set(['CZ2', 'CB', 'CG', 'CH2', 'CE3', 'CD1', 'CD2', 'CZ3', 'NE1', 'CE2']),
        'TYR': set(['CE1', 'OH', 'CB', 'CG', 'CZ', 'CD1', 'CD2', 'CE2']),
        'VAL': set(['CB', 'CG1', 'CG2']),
    }

    namedsele = cmd.get_unused_name('_')
    cmd.select(namedsele, selection, 0)

    namelists = defaultdict(list)
    cmd.iterate('(%s) and polymer' % namedsele,
            'namelists[model,segi,chain,resn,resi,resv].append(name)',
            space=locals())

    sele_dict = defaultdict(list)
    tmp_name = cmd.get_unused_name('_')

    for key, namelist in namelists.items():
        resn = key[3]
        if resn not in reference:
            if not quiet:
                print(' Unknown residue: + ' + resn)
            continue
        if not reference[resn].issubset(namelist):
            try:
                frag = fragments.get(resn.lower())
                for a in frag.atom:
                    a.segi = key[1]
                    a.chain = key[2]
                    a.resi = key[4]
                    a.resi_number = key[5]
                cmd.load_model(frag, tmp_name, 1, zoom=0)

                skey = '/%s/%s/%s/%s`%s' % key[:5]
                cmd.remove(skey + ' and not name N+C+O+OXT+CA')
                cmd.align(tmp_name, skey + ' and name N+C+CA', cycles=0)
                cmd.remove(tmp_name + ' and (name N+C+O+CA or hydro)')
                cmd.fuse('name CB and ' + tmp_name, 'name CA and ' + skey, move=0)
                if resn == 'PRO':
                    cmd.bond(skey + '/N', skey + '/CD')
                cmd.unpick()
                cmd.delete(tmp_name)

                sele_dict[key[0]].append(skey)

                if not quiet:
                    print(' Mutated ' + skey)
            except:
                print(' Mutating ' + skey + ' failed')

    for model in sele_dict:
        cmd.sort(model)
        sculpt_relax(' '.join(sele_dict[model]), 0, 0, model, cycles)

    cmd.delete(namedsele)

def peptide_rebuild(name, selection='all', cycles=1000, state=1, quiet=1):
    '''
DESCRIPTION

    Rebuild the peptide from selection. All atoms which are present in
    selection will be kept fixed, while atoms missing in selection are
    placed by sculpting.

USAGE

    peptide_rebuild name [, selection [, cycles [, state ]]]

SEE ALSO

    stub2ala, add_missing_atoms, peptide_rebuild_modeller
    '''
    from chempy import fragments, feedback, models

    cycles, state, quiet = int(cycles), int(state), int(quiet)

    # suppress feedback for model merging
    feedback['actions'] = False

    # work with named selection
    namedsele = cmd.get_unused_name('_')
    cmd.select(namedsele, '{} & present'.format(selection), 0)

    identifiers = []
    cmd.iterate(namedsele + ' and polymer and guide and alt +A',
            'identifiers.append([segi,chain,resi,resv,resn])', space=locals())

    model = models.Indexed()
    for (segi,chain,resi,resv,resn) in identifiers:
        try:
            fname = resn.lower() if resn != 'MSE' else 'met'
            frag = fragments.get(fname)
        except IOError:
            print(' Warning: unknown residue: ' + resn)
            continue

        for a in frag.atom:
            a.segi = segi
            a.chain = chain
            a.resi = resi
            a.resi_number = resv
            a.resn = resn

        model.merge(frag)

    if not quiet:
        print(' Loading model...')

    cmd.load_model(model, name, 1, zoom=0)
    if cmd.get_setting_boolean('auto_remove_hydrogens'):
        cmd.remove(name + ' and hydro')

    cmd.protect(name + ' in ' + namedsele)
    cmd.sculpt_activate(name)
    cmd.update(name, namedsele, 1, state)
    cmd.delete(namedsele)

    if not quiet:
        print(' Sculpting...')

    cmd.set('sculpt_field_mask', 0x003, name) # bonds and angles only
    cmd.sculpt_iterate(name, 1, int(cycles / 4))

    cmd.set('sculpt_field_mask', 0x09F, name) # local + torsions
    cmd.sculpt_iterate(name, 1, int(cycles / 4))

    cmd.set('sculpt_field_mask', 0x0FF, name) # ... + vdw
    cmd.sculpt_iterate(name, 1, int(cycles / 2))

    cmd.sculpt_deactivate(name)
    cmd.deprotect(name)
    cmd.unset('sculpt_field_mask', name)

    if not quiet:
        print(' Connecting peptide...')

    pairs = cmd.find_pairs(name + ' and name C', name + ' and name N', 1, 1, 2.0)
    for pair in pairs:
        cmd.bond(*pair)
    cmd.h_fix(name)

    if not quiet:
        print(' peptide_rebuild: done')


def get_seq(selection, chainbreak='/', unknown='A'):
    '''Gets the one-letter sequence, including residues without coordinates
    '''
    seq_list = []

    cmd.iterate('(%s) & polymer' % (selection),
            'seq_list.append((resn, resv))',
            space={'seq_list': seq_list})

    def seqbuilder():
        from . import one_letter
        prev_resv = None
        for resn, resv in seq_list:
            if resv != prev_resv:
                if prev_resv is not None and resv != prev_resv + 1:
                    yield chainbreak
                if resn in one_letter:
                    yield one_letter[resn]
                else:
                    print('Warning: unknown residue "%s"' % (resn))
                    yield unknown
                prev_resv = resv

    return ''.join(seqbuilder())


def peptide_rebuild_modeller(name, selection='all', hetatm=0, sequence=None,
        nmodels=1, hydro=0, quiet=1):
    '''
DESCRIPTION

    Remodel the given selection using modeller. This is useful for example to
    build incomplete sidechains. More complicated modelling tasks are not
    the intention of this simple interface.

    Side effects: Alters "type" property for MSE residues in selection
    (workaround for bug #3512313).

USAGE

    peptide_rebuild_modeller name [, selection [, hetatm [, sequence ]]]

ARGUMENTS

    name = string: new object name

    selection = string: atom selection

    hetatm = 0/1: read and model HETATMs (ligands) {default: 0}

    sequence = string: if provided, use this sequence instead of the
    template sequence {default: None}

    nmodels = int: number of models (states) to generate {default: 1}
    '''
    try:
        import modeller
        from modeller.automodel import automodel, allhmodel
    except ImportError:
        print(' Error: failed to import "modeller"')
        raise CmdException

    import tempfile, shutil, os
    from .editing import update_identifiers

    nmodels, hetatm, quiet = int(nmodels), int(hetatm), int(quiet)

    if int(hydro):
        automodel = allhmodel

    tempdir = tempfile.mkdtemp()
    pdbfile = os.path.join(tempdir, 'template.pdb')
    alnfile = os.path.join(tempdir, 'aln.pir')

    cwd = os.getcwd()
    os.chdir(tempdir)

    if not quiet:
        print(' Notice: PWD=%s' % (tempdir))

    try:
        modeller.log.none()
        env = modeller.environ()
        env.io.hetatm = hetatm

        # prevent PyMOL to put TER records before MSE residues (bug #3512313)
        cmd.alter('(%s) and polymer' % (selection), 'type="ATOM"')

        cmd.save(pdbfile, selection)
        mdl = modeller.model(env, file=pdbfile)

        aln = modeller.alignment(env)
        aln.append_model(mdl, align_codes='foo', atom_files=pdbfile)

        # get sequence from non-present atoms
        if not sequence and cmd.count_atoms('(%s) & !present' % (selection)):
            sequence = get_seq(selection)

        if sequence:
            aln.append_sequence(sequence)
            aln[-1].code = 'bar'
            aln.malign()
        aln.write(alnfile)

        a = automodel(env, alnfile=alnfile, sequence=aln[-1].code,
                knowns=[s.code for s in aln if s.prottyp.startswith('structure')])
        a.max_ca_ca_distance = 30.0

        if nmodels > 1:
            a.ending_model = nmodels
            from multiprocessing import cpu_count
            ncpu = min(cpu_count(), nmodels)
            if ncpu > 1:
                from modeller import parallel
                job = parallel.job(parallel.local_slave()
                        for _ in range(ncpu))
                a.use_parallel_job(job)

        a.make()

        for output in a.outputs:
            cmd.load(output['name'], name, quiet=quiet)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tempdir)

    cmd.align(name, selection, cycles=0)
    if not sequence:
        update_identifiers(name, selection)

    if not quiet:
        print(' peptide_rebuild_modeller: done')

cmd.extend('mutate', mutate)
cmd.extend('mutate_all', mutate_all)
cmd.extend('sculpt_relax', sculpt_relax)
cmd.extend('add_missing_atoms', add_missing_atoms)
cmd.extend('peptide_rebuild', peptide_rebuild)
cmd.extend('peptide_rebuild_modeller', peptide_rebuild_modeller)

cmd.auto_arg[0].update([
    ('mutate', cmd.auto_arg[0]['align']),
    ('mutate_all', cmd.auto_arg[0]['align']),
    ('sculpt_relax', cmd.auto_arg[0]['align']),
    ('add_missing_atoms', cmd.auto_arg[0]['zoom']),
])
cmd.auto_arg[1].update([
    ('peptide_rebuild', cmd.auto_arg[0]['zoom']),
    ('peptide_rebuild_modeller', cmd.auto_arg[0]['zoom']),
])

# vi: expandtab:smarttab
