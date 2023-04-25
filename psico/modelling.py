'''
(c) 2010-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

_auto_arg0_align = cmd.auto_arg[0]['align']
_auto_arg1_align = cmd.auto_arg[1]['align']


def _assert_package_import():
    if not __name__.endswith('.modelling'):
        raise CmdException("Must do 'import psico.modelling' instead of 'run ...'")

def mutate(selection, new_resn, inplace=0, sculpt=0, hydrogens='auto', mode=3,
        quiet=1, *, _self=cmd):
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

    mode = 0: select rotamer with best clash score
    mode = 1: take chi angles from original residue
    mode = 2: first rotamer
    mode = 3: wizard default {default}

EXAMPLE

    fetch 2x19, async=0
    select x, A/CYS`122/
    zoom x
    mutate x, LYS
    '''
    from pymol.wizard import mutagenesis
    _assert_package_import()
    from . import three_letter
    from . import selecting

    inplace, sculpt = int(inplace), int(sculpt)
    mode = int(mode)
    quiet = int(quiet)
    org = _self.get_object_list(selection)[0]
    tmp = _self.get_unused_name()
    new_resn = new_resn.upper()
    new_resn = three_letter.get(new_resn, new_resn)

    if inplace:
        cpy = org
    else:
        cpy = _self.get_unused_name(org + '_cpy')
        _self.create(cpy, org, -1, 1, zoom=0)

    scr = []
    _self.iterate('first (%s)' % selection, 'scr[:] = (segi,chain,resi,resn)', space={'scr': scr})
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
                old_chi.append(_self.get_dihedral(*args))
            except:
                break

    _self.remove('%s and not name CA+C+N+O+OXT' % (res))

    # start the wizard to count the number of rotamers for this residue
    _self.wizard("mutagenesis")
    _self.get_wizard().set_mode(new_resn)
    _self.get_wizard().set_hyd(hydrogens)

    with selecting.select_temporary(res, _self=_self) as res_named_sele:
        _self.get_wizard().do_select(res_named_sele)

    def get_best_state_bump():
        best_state = (1, 1e9)
        _self.create(tmp, '%s and not name CA+C+N+O or (%s within 8.0 of (%s and name CB))' % \
                (mutagenesis.obj_name, cpy, mutagenesis.obj_name), zoom=0, singletons=1)
        _self.bond('name CB and %s in %s' % (tmp, mutagenesis.obj_name),
                'name CA and %s in %s' % (tmp, res))
        _self.sculpt_activate(tmp)
        for i in range(1, _self.count_states(tmp)+1):
            score = _self.sculpt_iterate(tmp, state=i)
            if not quiet:
                print('Frame %d Score %.2f' % (i, score))
            if score < best_state[1]:
                best_state = (i, score)
        _self.delete(tmp)
        if not quiet:
            print(' Best: Frame %d Score %.2f' % best_state)
        return best_state

    if _self.count_states(mutagenesis.obj_name) < 2 or mode > 0:
        best_state = (1, -1.0)
    else:
        best_state = get_best_state_bump()

    if mode != 3:
        _self.frame(best_state[0])

    _self.get_wizard().apply()
    _self.wizard()

    if mode == 1:
        atoms = [res + '/' + name for name in chi_atoms.get(new_resn, [])]
        for args in zip(atoms, atoms[1:], atoms[2:], atoms[3:], old_chi):
            _self.set_dihedral(*args)
        _self.unpick()

    if sculpt > 0:
        sculpt_relax(res, 0, sculpt == 2, cpy)

    return cpy

def mutate_all(selection, new_resn, inplace=1, sculpt=0, *args, _self=cmd, **kwargs):
    '''
DESCRIPTION

    Mutate all residues in selection. By default do mutation in-place (unlike
    the 'mutate' command which by default works on a copy).

    FOR SCULPTING ONLY SUPPORTS SELECTIONS WITHIN THE SAME MODEL!

SEE ALSO

    mutate
    '''
    inplace, sculpt = int(inplace), int(sculpt)

    if sculpt and len(_self.get_object_list('(' + selection + ')')) > 1:
        raise CmdException('Sculpting in multiple models not supported')

    kwargs.pop('_self', None)
    sele_list = set()
    _self.iterate(selection,
            'sele_list.add("/%s/%s/%s/%s" % (model, segi, chain, resi))',
            space={'sele_list': sele_list})
    for sele in sele_list:
        mutate(sele, new_resn, inplace, sculpt and not inplace, *args, **kwargs)
    if sculpt and inplace:
        sculpt_relax('(' + ' '.join(sele_list) + ')', 0, sculpt == 2)

def sculpt_relax(selection, backbone=1, neighbors=0, model=None, cycles=100,
        state=0, quiet=1, *, _self=cmd):
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
    org = _self.get_object_list(sele)[0]
    if model is None:
        model = org
    elif model != org:
        sele = sele.replace('(%s)' % org, '(%s)' % model)

    _self.protect()
    _self.deprotect(sele)
    if not backbone:
        _self.protect('name CA+C+N+O+OXT')

    _self.sculpt_activate(model, state)
    _self.set('sculpt_vdw_weight', 0.25, model) # Low VDW forces
    _self.set('sculpt_field_mask', 0x1FF, model) # Default

    if neighbors:
        _self.sculpt_iterate(model, state, int(cycles * 0.25))
        _self.deprotect('byres (%s within 6.0 of (%s))' % (model, sele))
        if not backbone:
            _self.protect('name CA+C+N+O+OXT')
        _self.sculpt_iterate(model, state, cycles=int(cycles * 0.50))
    else:
        _self.sculpt_iterate(model, state, int(cycles * 0.75))

    _self.set('sculpt_field_mask', 0x01F, model) # Local Geometry Only
    _self.sculpt_iterate(model, state, int(cycles * 0.25))

    _self.unset('sculpt_vdw_weight', model)
    _self.unset('sculpt_field_mask', model)
    _self.sculpt_deactivate(model)
    _self.deprotect()

def add_missing_atoms(selection='all', cycles=200, quiet=1, *, _self=cmd):
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

    namedsele = _self.get_unused_name('_')
    _self.select(namedsele, selection, 0)

    namelists = defaultdict(list)
    _self.iterate('(%s) and polymer' % namedsele,
            'namelists[model,segi,chain,resn,resi,resv].append(name)',
            space=locals())

    sele_dict = defaultdict(list)
    tmp_name = _self.get_unused_name('_')

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
                _self.load_model(frag, tmp_name, 1, zoom=0)

                skey = '/%s/%s/%s/%s`%s' % key[:5]
                _self.remove(skey + ' and not name N+C+O+OXT+CA')
                _self.align(tmp_name, skey + ' and name N+C+CA', cycles=0)
                _self.remove(tmp_name + ' and (name N+C+O+CA or hydro)')
                _self.fuse('name CB and ' + tmp_name, 'name CA and ' + skey, move=0)
                if resn == 'PRO':
                    _self.bond(skey + '/N', skey + '/CD')
                _self.unpick()
                _self.delete(tmp_name)

                sele_dict[key[0]].append(skey)

                if not quiet:
                    print(' Mutated ' + skey)
            except:
                print(' Mutating ' + skey + ' failed')

    for model in sele_dict:
        _self.sort(model)
        sculpt_relax(' '.join(sele_dict[model]), 0, 0, model, cycles)

    _self.delete(namedsele)

def peptide_rebuild(name, selection='all', cycles=1000, state=1, quiet=1, *, _self=cmd):
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
    namedsele = _self.get_unused_name('_')
    _self.select(namedsele, '{} & present'.format(selection), 0)

    identifiers = []
    _self.iterate(namedsele + ' and polymer and guide and alt +A',
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

    _self.load_model(model, name, 1, zoom=0)
    if _self.get_setting_boolean('auto_remove_hydrogens'):
        _self.remove(name + ' and hydro')

    _self.protect(name + ' in ' + namedsele)
    _self.sculpt_activate(name)
    _self.update(name, namedsele, 1, state)
    _self.delete(namedsele)

    if not quiet:
        print(' Sculpting...')

    _self.set('sculpt_field_mask', 0x003, name) # bonds and angles only
    _self.sculpt_iterate(name, 1, int(cycles / 4))

    _self.set('sculpt_field_mask', 0x09F, name) # local + torsions
    _self.sculpt_iterate(name, 1, int(cycles / 4))

    _self.set('sculpt_field_mask', 0x0FF, name) # ... + vdw
    _self.sculpt_iterate(name, 1, int(cycles / 2))

    _self.sculpt_deactivate(name)
    _self.deprotect(name)
    _self.unset('sculpt_field_mask', name)

    if not quiet:
        print(' Connecting peptide...')

    pairs = _self.find_pairs(name + ' and name C', name + ' and name N', 1, 1, 2.0)
    for pair in pairs:
        _self.bond(*pair)
    _self.h_fix(name)

    if not quiet:
        print(' peptide_rebuild: done')


def get_seq(selection, chainbreak='/', unknown='A', *, _self=cmd):
    '''Gets the one-letter sequence, including residues without coordinates
    '''
    seq_list = []

    _self.iterate('(%s) & polymer' % (selection),
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
        nmodels=1, hydro=0, quiet=1, *, _self=cmd):
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
    import modeller
    from modeller.automodel import automodel, allhmodel

    import tempfile, shutil, os
    _assert_package_import()
    from .editing import update_identifiers

    nmodels, hetatm, quiet = int(nmodels), int(hetatm), int(quiet)

    if int(hydro):
        automodel = allhmodel  # noqa: F811 Redefinition of unused

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
        _self.alter('(%s) and polymer' % (selection), 'type="ATOM"')

        _self.save(pdbfile, selection)
        mdl = modeller.model(env, file=pdbfile)

        aln = modeller.alignment(env)
        aln.append_model(mdl, align_codes='foo', atom_files=pdbfile)

        # get sequence from non-present atoms
        if not sequence and _self.count_atoms('(%s) & !present' % (selection)):
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
            _self.load(output['name'], name, quiet=quiet)
    finally:
        os.chdir(cwd)
        shutil.rmtree(tempdir)

    _self.align(name, selection, cycles=0)
    if not sequence:
        update_identifiers(name, selection, _self=_self)

    if not quiet:
        print(' peptide_rebuild_modeller: done')


@cmd.extendaa(_auto_arg0_align, _auto_arg1_align)
def update_align(mobile: str,
                 target: str,
                 state: int = 1,
                 *,
                 fix: str = "none",
                 quiet: int = 1,
                 _self=cmd):
    """
DESCRIPTION

    Update (and optionally fix) coordinates based on sequence alignment.
    """
    aln = _self.get_unused_name("aln_hom")
    _self.align(mobile, target, object=aln, cycles=0, max_gap=-1)
    try:
        mobile_aln = f"({mobile}) & {aln}"
        _self.update(mobile_aln,
                     f"({target}) & {aln}",
                     state,
                     state,
                     matchmaker=0,
                     quiet=quiet)
        if fix == "restrain":
            _self.reference("store", mobile_aln, state, quiet=quiet)
            _self.flag("restrain", mobile_aln, "set", quiet=quiet)
        elif fix == "fix":
            _self.flag("fix", mobile_aln, "set", quiet=quiet)
        elif fix == "protect":
            _self.protect(mobile_aln, quiet=quiet)
        elif fix != "none":
            raise ValueError(fix)
    finally:
        _self.delete(aln)


@cmd.extendaa(_auto_arg0_align, _auto_arg1_align)
def sculpt_homolog(mobile: str,
                   target: str,
                   state: int = 1,
                   cycles: int = 1000,
                   *,
                   fix: str = "restrain",
                   quiet: int = 1,
                   _self=cmd):
    """
DESCRIPTION

    Sculpt mobile towards target, based on sequence alignment.

ARGUMENTS

    cycles: Number of sculpt iterations

    fix = restrain | fix | protect | none: Method for fixing updated atoms
    """
    (mobile_object, ) = _self.get_object_list(mobile)
    _self.sculpt_activate(mobile_object, state)
    update_align(mobile, target, state, fix=fix, quiet=quiet, _self=_self)
    _self.sculpt_iterate(mobile, state, cycles)


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
