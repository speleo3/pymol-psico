import psico.editing
import os
import pytest
from pymol import cmd
from pytest import approx

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
FILENAME_MULTISTATE = os.path.join(DATA_DIR, '1nmr-frag-nohydro.pdb')
FILENAME_WITH_MSE = os.path.join(DATA_DIR, '2x19-frag-mse.pdb')


def test_split_chains():
    cmd.reinitialize()
    cmd.fab('ACD', 'm1', chain='A')
    cmd.fab('EFG', 'm2', chain='B')
    cmd.create('m3', 'm1 m2')
    psico.editing.split_chains('m3')
    assert cmd.get_chains('m3_A') == ['A']
    assert cmd.get_chains('m3_B') == ['B']
    psico.editing.split_chains('m3', 'foo_')
    assert cmd.get_chains('foo_0001') == ['A']


def test_split_molecules():
    cmd.reinitialize()
    cmd.fab('ACD', 'm1')
    cmd.fab('EFG', 'm2')
    cmd.create('m3', 'm1 m2')
    psico.editing.split_molecules('m3')
    assert cmd.count_atoms('mol_01') == cmd.count_atoms('m1')
    assert cmd.count_atoms('mol_02') == cmd.count_atoms('m2')
    psico.editing.split_molecules('m3', 'foo_')
    assert cmd.count_atoms('foo_01') == cmd.count_atoms('m1')


def test_split():
    cmd.reinitialize()
    cmd.fab('ACD', 'm1')
    psico.editing.split('byres', 'm1')
    assert cmd.count_atoms('entity01') == cmd.count_atoms('m1 & resn ALA')
    assert cmd.count_atoms('entity02') == cmd.count_atoms('m1 & resn CYS')
    assert cmd.count_atoms('entity03') == cmd.count_atoms('m1 & resn ASP')
    psico.editing.split('byres', 'm1', 'foo_')
    assert cmd.count_atoms('foo_01') == cmd.count_atoms('m1 & resn ALA')


def test_rmsf2b():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.editing.rmsf2b()
    my_values = []
    cmd.iterate('resi 80 & guide', 'my_values.append(b)', space=locals())
    assert my_values == approx([3.266], rel=1e-2)
    psico.editing.rmsf2b(linearscale=-1, var='q')
    my_values = []
    cmd.iterate('resi 80 & guide', 'my_values.append(q)', space=locals())
    assert my_values == approx([842.231], rel=1e-2)


def test_set_sequence():
    cmd.reinitialize()
    cmd.fab('ACD', 'm1')
    psico.editing.set_sequence('EFG')
    myseq = []
    cmd.iterate('guide', 'myseq.append(resn)', space=locals())
    assert myseq == ['GLU', 'PHE', 'GLY']
    psico.editing.set_sequence('HI', start=2)
    myseq = []
    cmd.iterate('guide', 'myseq.append(resn)', space=locals())
    assert myseq == ['GLU', 'HIS', 'ILE']


def test_alphatoall():
    cmd.reinitialize()
    cmd.fab('ACD', 'm1')
    cmd.alter('guide', 'color = resv')
    psico.editing.alphatoall('all', 'color')
    my_values = []
    cmd.iterate('all', 'my_values.append(color)', space=locals())
    assert set(my_values) == set([1, 2, 3])
    cmd.alter('name N', 'color = resv + 3')
    psico.editing.alphatoall('all', 'color', operator='(name N) and')
    my_values = []
    cmd.iterate('all', 'my_values.append(color)', space=locals())
    assert set(my_values) == set([4, 5, 6])


def test_mse2met():
    cmd.reinitialize()
    cmd.load(FILENAME_WITH_MSE)
    assert cmd.count_atoms('resn MSE') == 16
    assert cmd.count_atoms('name SE') == 2
    psico.editing.mse2met()
    assert cmd.count_atoms('resn MSE') == 0
    assert cmd.count_atoms('name SE') == 0
    assert cmd.count_atoms('resn MET') == 16
    assert cmd.count_atoms('name SD') == 2


def test_polyala():
    cmd.reinitialize()
    cmd.fab('EFG', 'm1')
    psico.editing.polyala()
    assert cmd.count_atoms('name CA') == 3
    assert cmd.count_atoms('name CB') == 2
    assert cmd.count_atoms('name CG+CD+CD1+CD2+2HG+3HG+1HD+2HD') == 0


def test_stub2ala():
    cmd.reinitialize()
    cmd.fab('EFG', 'm1')
    cmd.remove('not (backbone or name CB)')
    psico.editing.stub2ala()
    my_values = []
    cmd.iterate('guide', 'my_values.append(resn)', space=locals())
    assert my_values == ['ALA', 'ALA', 'GLY']


def test_remove_alt():
    cmd.reinitialize()
    cmd.fab('AC', 'm1')
    cmd.alter('resn ALA', 'alt="A"')
    cmd.alter('resn CYS', 'alt="B"')
    cmd.create('m2', 'm1')
    psico.editing.remove_alt('m1')
    psico.editing.remove_alt('m2', keep='B')
    assert cmd.count_atoms('m1') == 10
    assert cmd.count_atoms('m2') == 11


@pytest.mark.exe
def test_dssp():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.editing.dssp(raw='resn')
    my_values = []
    cmd.iterate('guide', 'my_values.append((ss, resn))', space=locals())
    assert my_values == [\
        ('L', ' '), ('L', 'T'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('L', 'T'), ('L', 'T'),
        ('L', ' '), ('L', ' '), ('L', 'T'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('L', 'T'), ('L', 'T'), ('L', ' '),
        ('L', ' ')]


@pytest.mark.exe
def test_stride():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.editing.stride(raw='resn')
    my_values = []
    cmd.iterate('guide', 'my_values.append((ss, resn))', space=locals())
    assert my_values == [\
        ('L', 'C'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('L', 'C'),
        ('L', 'C'), ('L', 'C'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('L', 'T'), ('L', 'T'), ('L', 'T'),
        ('L', 'C')]


@pytest.mark.web
def test_sst():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.editing.sst(raw='resn')
    my_values = []
    cmd.iterate('guide', 'my_values.append((ss, resn))', space=locals())
    assert my_values == [\
        ('L', 'C'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'G'), ('H', 'G'), ('H', 'G'), ('L', 'C'),
        ('L', 'C'), ('H', 'I'), ('H', 'I'), ('H', 'I'), ('H', 'I'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'), ('H', 'H'),
        ('H', 'H'), ('L', 'C'), ('L', '3'), ('L', '3'), ('L', '3'),
        ('L', '3')]


def test_set_phipsi():
    cmd.reinitialize()
    cmd.fab('ACDEF', 'm1')
    psico.editing.set_phipsi('all', 160, 160)
    psico.editing.set_phipsi('resi 2+4', 120)
    psico.editing.set_phipsi('resi 3+4', psi=-140)
    assert cmd.get_phipsi('resi 2')['m1', 12] == approx((120., 160.))
    assert cmd.get_phipsi('resi 3')['m1', 23] == approx((160., -140.))
    assert cmd.get_phipsi('resi 4')['m1', 35] == approx((120., -140.))


def test_update_identifiers():
    cmd.reinitialize()
    cmd.fab('ACD', 'm1')
    cmd.fab('ACD', 'm2')
    cmd.remove('m1 and not backbone')
    cmd.alter('m1', '(segi, chain) = ("Segi", "Chain")')
    psico.editing.update_identifiers('m2', 'm1', identifiers='segi chain')
    assert cmd.get_chains('m2') == ['Chain']
    my_values = []
    cmd.iterate('m2', 'my_values.append(segi)', space=locals())
    assert set(my_values) == set(["Segi"])
