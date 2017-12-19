'''
(c) Thomas Holder, Schrodinger Inc.
'''

from pymol import cmd, CmdException


def get_fixed_indices(selection, state, _self):
    fixed_list = []
    _self.iterate_state(state, selection,
            '_append(flags & 0x8)',
            space={'_append': fixed_list.append})
    return [idx for (idx, fixed) in enumerate(fixed_list) if fixed]


def load_or_update(molstr, name, sele, state, _self):
    update = not name

    if update:
        name = _self.get_unused_name('_minimized')
    else:
        _self.delete(name)

    _self.load(molstr, name, 1, 'molstr', zoom=0)

    try:
        from psico.fitting import xfit
        xfit(name, sele, 1, state, match='none', cycles=100, guide=0)
    except Exception as e:
        print('xfit failed, fallback to cmd.fit')
        _self.fit(name, sele, 1, state, cycles=5, matchmaker=-1)

    if update:
        _self.update(sele, name, state, 1, matchmaker=0)
        _self.delete(name)


def minimize_ob(selection='enabled', state=-1, ff='UFF', nsteps=500,
        conv=0.0001, cutoff=0, cut_vdw=6.0, cut_elec=8.0,
        name='', quiet=1, _self=cmd):
    '''
DESCRIPTION

    Emergy minimization with openbabel

    Supports fixed atoms (flag fix)

ARGUMENTS

    selection = str: atom selection

    state = int: object state {default: -1}

    ff = GAFF|MMFF94s|MMFF94|UFF|Ghemical: force field {default: UFF}

    nsteps = int: number of steps {default: 500}
    '''
    import openbabel as ob

    state = int(state)

    sele = _self.get_unused_name('_sele')
    _self.select(sele, selection, 0)

    try:
        ioformat = 'mol'
        molstr = _self.get_str(ioformat, sele, state)

        obconversion = ob.OBConversion()
        obconversion.SetInAndOutFormats(ioformat, ioformat)

        mol = ob.OBMol()
        obconversion.ReadString(mol, molstr)

        # add hydrogens
        orig_ids = [a.GetId() for a in ob.OBMolAtomIter(mol)]
        mol.AddHydrogens()
        added_ids = set(a.GetId() for a in ob.OBMolAtomIter(mol)).difference(orig_ids)

        consttrains = ob.OBFFConstraints()
        consttrains.Setup(mol)

        # atoms with "flag fix"
        fixed_indices = get_fixed_indices(sele, state, _self)
        for idx in fixed_indices:
            consttrains.AddAtomConstraint(idx + 1)

        # setup forcefield (one of: GAFF, MMFF94s, MMFF94, UFF, Ghemical)
        ff = ob.OBForceField.FindForceField(ff)
        ff.Setup(mol, consttrains)

        if int(cutoff):
            ff.EnableCutOff(True)
            ff.SetVDWCutOff(float(cut_vdw))
            ff.SetElectrostaticCutOff(float(cut_elec))

        # run minimization
        ff.SteepestDescent(int(nsteps) // 2, float(conv))
        ff.ConjugateGradients(int(nsteps) // 2, float(conv))
        ff.GetCoordinates(mol)

        # remove previously added hydrogens
        for hydro_id in added_ids:
            mol.DeleteAtom(mol.GetAtomById(hydro_id))

        molstr = obconversion.WriteString(mol)
        load_or_update(molstr, name, sele, state, _self)

        if not int(quiet):
            print(' Energy: %8.2f %s' % (ff.Energy(), ff.GetUnit()))
    finally:
        _self.delete(sele)


def minimize_rdkit(selection='enabled', state=-1, ff='MMFF94', nsteps=200,
        name='', quiet=1, _self=cmd):
    '''
DESCRIPTION

    Emergy minimization with RDKit

    Supports fixed atoms (flag fix)

ARGUMENTS

    selection = str: atom selection

    state = int: object state {default: -1}

    ff = MMFF94s|MMFF94|UFF: force field {default: MMFF94}

    nsteps = int: number of steps {default: 200}
    '''
    from rdkit import Chem
    from rdkit.Chem import AllChem

    state = int(state)

    sele = _self.get_unused_name('_sele')
    _self.select(sele, selection, 0)

    try:
        molstr = _self.get_str('mol', sele, state)
        mol = Chem.MolFromMolBlock(molstr, True, False)

        if mol is None:
            raise CmdException('Failed to load molecule into RDKit. '
                    'Please check bond orders and formal charges.')

        # setup forcefield
        if ff.startswith('MMFF'):
            ff = AllChem.MMFFGetMoleculeForceField(mol,
                    AllChem.MMFFGetMoleculeProperties(mol, ff,
                        0 if int(quiet) else 1))
        elif ff == 'UFF':
            ff = AllChem.UFFGetMoleculeForceField(mol)
        else:
            raise CmdException('unknown forcefield: ' + ff)

        if ff is None:
            raise CmdException('forcefield setup failed')

        # atoms with "flag fix"
        for idx in get_fixed_indices(sele, state, _self):
            ff.AddFixedPoint(idx)

        # run minimization
        if ff.Minimize(int(nsteps)) != 0:
            print(" Warning: minimization did not converge")

        molstr = Chem.MolToMolBlock(mol)
        load_or_update(molstr, name, sele, state, _self)

        if not int(quiet):
            print(' Energy: %8.2f %s' % (ff.CalcEnergy(), 'kcal/mol'))
    finally:
        _self.delete(sele)


cmd.extend('minimize_ob', minimize_ob)
cmd.extend('minimize_rdkit', minimize_rdkit)

cmd.auto_arg[0].update([
    ('minimize_ob', cmd.auto_arg[0]['zoom']),
    ('minimize_rdkit', cmd.auto_arg[0]['zoom']),
])

# vi:expandtab:smarttab
