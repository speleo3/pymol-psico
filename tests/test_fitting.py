import psico.fitting
import psico.querying
import os
from pymol import cmd
import pytest
from pytest import approx

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
FILENAME_MULTISTATE = os.path.join(DATA_DIR, '1nmr-frag-nohydro.pdb')

skipif_GIL_bug = pytest.mark.skipif(
    psico.pymol_version_tuple < (2, 5, 6),
    reason="Affected by transform_object GIL bug, "
    "fixed in https://github.com/schrodinger/pymol-open-source/pull/285")


def test_alignwithanymethod():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    cmd.create("m3", "m1", 3, 1)
    psico.fitting.alignwithanymethod("m3", "m1", "align super", target_state=1, async_=0)
    assert cmd.get_object_list() == ["m1", "m3", "m3_align01", "m3_super01"]


@skipif_GIL_bug
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


def _assert_object_matrix_equal(obj, state, matrix):
    assert cmd.get_object_matrix(obj, state, 0) == approx(matrix, abs=1e-3)


@skipif_GIL_bug
@pytest.mark.exe
def test_theseus():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    cmd.create("m3", "m1", 3, 1)
    psico.fitting.theseus("m3", "m1")
    _assert_object_matrix_equal("m3", 1, [
        0.9992, 0.0347, -0.0222, -0.1863, -0.0349, 0.9993, -0.01, -0.2639,
        0.0218, 0.0108, 0.9997, 0.3903, 0.0, 0.0, 0.0, 1.0
    ])


@skipif_GIL_bug
@pytest.mark.exe
def test_intra_theseus():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE, "m1")
    psico.fitting.intra_theseus("m1")
    _assert_object_matrix_equal("m1", 1, [
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 1.0
    ])
    _assert_object_matrix_equal("m1", 2, [
        0.99876, 0.02313, -0.04404, -0.09647, -0.02267, 0.99968, 0.01101,
        -0.34898, 0.04428, -0.01, 0.99897, 0.61188, 0.0, 0.0, 0.0, 1.0
    ])
    _assert_object_matrix_equal("m1", 3, [
        0.99968, 0.02362, -0.00903, -0.10755, -0.02368, 0.99969, -0.0074,
        -0.16009, 0.00886, 0.00762, 0.99993, 0.29786, 0.0, 0.0, 0.0, 1.0
    ])
    _assert_object_matrix_equal("m1", 4, [
        0.99656, 0.08201, -0.01205, -0.21902, -0.08149, 0.99592, 0.03863,
        -0.96714, 0.01517, -0.03752, 0.99918, 0.50427, 0.0, 0.0, 0.0, 1.0
    ])
    psico.fitting.intra_theseus("m1", cov=1, cycles=50)
    _assert_object_matrix_equal("m1", 1, [
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 1.0
    ])
    _assert_object_matrix_equal("m1", 2, [
        0.99957, 0.01668, -0.02415, -0.02072, -0.01597, 0.99945, 0.02918,
        -0.26491, 0.02463, -0.02878, 0.99928, 0.5978, 0.0, 0.0, 0.0, 1.0
    ])
    _assert_object_matrix_equal("m1", 3, [
        0.99933, 0.01346, 0.03418, -0.04087, -0.01223, 0.99928, -0.03588,
        -0.0714, -0.03463, 0.03544, 0.99877, -0.34112, 0.0, 0.0, 0.0, 1.0
    ])
    _assert_object_matrix_equal("m1", 4, [
        0.99649, 0.07346, 0.04025, -0.15404, -0.07444, 0.99695, 0.02331,
        -0.88119, -0.03841, -0.02622, 0.99892, -0.13483, 0.0, 0.0, 0.0, 1.0
    ])


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


# def dyndom_parse_info(filename, selection='(all)', quiet=0):
# def dyndom(mobile, target, window=5, domain=20, ratio=1.0, exe='', transform=1,
# def matchmaker(mobile, target, match):
# def prosmart(mobile, target, mobile_state=1, target_state=1,
# def xfit(mobile, target, mobile_state=-1, target_state=-1, load_b=0,
# def promix(mobile, target, K=0, prefix=None, mobile_state=-1, target_state=-1,
# def intra_promix(selection, K=0, prefix=None, conformers=0, guide=1,
# def intra_boxfit(selection="polymer", center=[0.5, 0.5, 0.5], _self=cmd):
