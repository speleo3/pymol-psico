import psico.fitting
import psico.querying
import os
from pymol import cmd
import pytest
from pytest import approx

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
FILENAME_MULTISTATE = os.path.join(DATA_DIR, '1nmr-frag-nohydro.pdb')


@pytest.mark.skipif(
    psico.pymol_version_tuple < (2, 6),
    reason="Affected by transform_object GIL bug, "
    "fixed in https://github.com/schrodinger/pymol-open-source/pull/285")
@pytest.mark.exe
def test_tmalign():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    cmd.create("m2", "m1", 1, 1)
    score = psico.fitting.tmalign("m1",
                                  "m2",
                                  mobile_state=1,
                                  target_state=1,
                                  object="aln")
    assert score == approx(1.0)
    assert cmd.count_atoms("aln") == 62
    score = psico.fitting.tmalign("m1 & resi 51-72",
                                  "m2 & resi 62-82",
                                  mobile_state=3,
                                  target_state=1,
                                  object="aln")
    assert score == approx(0.50016)
    assert cmd.count_atoms("aln") == 22


def test_gdt_ts():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    cmd.create("m2", "m1", 1, 1)
    cmd.create("m3", "m1", 3, 1)
    r = psico.fitting.gdt_ts("m2", "m3")
    assert r == approx(0.74193548)


def test_get_rmsd_func(monkeypatch):
    X = [[0, 0, 0], [4, 0, 0], [0, 4, 0]]
    Y = [[0, 0, 0], [4, 0, 0], [0, 8, 0]]
    with monkeypatch.context() as monkeycontext:
        monkeycontext.delattr("chempy.cpv.fit")  # use numpy
        func = psico.fitting.get_rmsd_func()
        assert func(func.array(X), func.array(Y)) == approx(1.806, abs=0.001)
    with monkeypatch.context() as monkeycontext:
        monkeycontext.delattr("numpy.linalg.svd")  # use chempy
        func = psico.fitting.get_rmsd_func()
        assert func(func.array(X), func.array(Y)) == approx(1.806, abs=0.1)


def test_local_rms():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    r = psico.fitting.local_rms("m1", "m1", 10, 1, 3, match="none")
    ref = {
        52: 0.083, 53: 0.080, 54: 0.088, 55: 0.098, 56: 0.122, 57: 0.160,
        58: 0.181, 59: 0.175, 60: 0.164, 61: 0.164, 62: 0.160, 63: 0.140,
        64: 0.132, 65: 0.111, 66: 0.116, 67: 0.111, 68: 0.095, 69: 0.095,
        70: 0.160, 71: 0.223, 72: 0.243, 73: 0.250, 74: 0.566, 75: 1.669,
        76: 1.705, 77: 1.749, 78: 1.781, 79: 1.760, 80: 1.708, 81: 1.552,
        82: 1.668
    }
    assert r == approx(ref, abs=0.001)


def test_extra_fit():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    cmd.split_states("m1")
    psico.fitting.extra_fit("m1_*", object="aln")
    assert cmd.count_atoms("aln & guide") == 116
    assert cmd.count_atoms("aln & m1") == 0


def test_intra_xfit():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.fitting.intra_xfit("all", load_b=1, cycles=10)
    b_list = psico.querying.iterate_to_list("guide", "b")
    b_list_ref = [
        -3.932, -3.502, -3.806, -3.425, -4.862, -5.246, -4.358, -3.846, -4.348,
        -5.397, -4.353, -4.296, -4.037, -4.310, -4.874, -5.439, -4.552, -3.804,
        -3.987, -3.468, -3.019, -3.061, -1.961, -1.850, -1.721, -1.445, -1.088,
        -0.109, 2.399, 1.314, 2.389
    ]
    assert b_list == approx(b_list_ref, abs=1e-3)


def test_intra_center():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.fitting.intra_center("resi 80")
    assert cmd.intra_rms_cur("*") == approx(
        [-1.0, 1.425962, 8.287644, 7.676795])


# def alignwithanymethod(mobile, target, methods=None, async_=1, quiet=1, **kwargs):
# def dyndom_parse_info(filename, selection='(all)', quiet=0):
# def dyndom(mobile, target, window=5, domain=20, ratio=1.0, exe='', transform=1,
# def matchmaker(mobile, target, match):
# def theseus(mobile, target, match='align', cov=0, cycles=200,
# def intra_theseus(selection, state=1, cov=0, cycles=200,
# def prosmart(mobile, target, mobile_state=1, target_state=1,
# def xfit(mobile, target, mobile_state=-1, target_state=-1, load_b=0,
# def promix(mobile, target, K=0, prefix=None, mobile_state=-1, target_state=-1,
# def intra_promix(selection, K=0, prefix=None, conformers=0, guide=1,
# def intra_boxfit(selection="polymer", center=[0.5, 0.5, 0.5], _self=cmd):
