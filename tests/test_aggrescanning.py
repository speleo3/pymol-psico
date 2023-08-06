import psico.aggrescanning
import psico.querying
from pymol import cmd
from pytest import approx
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / 'data'


def approx3(data):
    return approx(data, abs=1e-3)


def test_aggrescan1d_scores():
    # peptidB from web example (truncated)
    p2 = psico.aggrescanning.aggrescan1d_scores("DAEFRHDSGYEVHHQKLVFFA")
    table = p2.pop("table")
    assert p2["a3vSA"] == approx(-0.139, abs=1e-3)
    assert p2["nHS"] == 1
    assert p2["NnHS"] == approx(4.762, abs=1e-3)
    assert p2["AAT"] == approx(4.727, abs=1e-3)
    assert p2["THSA"] == approx(4.605, abs=1e-3)
    assert p2["TA"] == approx(-2.601, abs=1e-3)
    assert p2["AATr"] == approx(0.225, abs=1e-3)
    assert p2["THSAr"] == approx(0.219, abs=1e-3)
    assert p2["Na4vSS"] == approx(-13.9, abs=1e-2)
    assert table["a4v"] == approx3([
        -0.631, -0.631, -0.554, -0.393, -0.753, -0.530, -0.988, -0.508, -0.584,
        0.102, -0.045, -0.145, -0.623, -0.527, -0.570, -0.044, 0.513, 1.110,
        1.289, 0.796, 0.796
    ])
    assert table["HSA"] == approx3([
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0.122, 0, 0, 0, 0, 0, 0, 4.605, 4.605,
        4.605, 4.605, 4.605
    ])
    assert table["NHSA"] == approx3([
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.921, 0.921, 0.921,
        0.921, 0.921
    ])
    assert table["a4vAHS"] == approx3([
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0.102, 0, 0, 0, 0, 0, 0, 0.901, 0.901,
        0.901, 0.901, 0.901
    ])


def test_aggrescan1d():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    psico.aggrescanning.aggrescan1d()
    b_list = psico.querying.iterate_to_list("guide", "b")
    b_list_ref = [
        0.3277, 0.3277, 0.3277, 0.537, 0.6041, 0.7573, 0.4206, 0.037, -0.2134,
        0.1799, 0.1799, 0.1799, 0.0009, 0.305, 0.2363, 0.2113, -0.091, -0.3306,
        -0.5504, -0.0884, -0.5179, -0.5021, -0.2267, -0.3177, -0.0606, -0.1709,
        -0.1709, -0.145, -0.2213, -0.582, -0.6507, -0.722, -0.3287, -0.6367,
        -0.5086, -0.5086, -0.5514, -0.5256, -0.6263, -0.6894, -0.3814, -0.0831,
        0.174, 0.2734, 0.3163, 0.3176, 0.3176, -0.1444, -0.6573, -0.7286,
        -0.8293, -0.719, -0.346, -0.5851, -0.6457, -0.2179, -0.3274, 0.11,
        -0.0431, -0.3733, -0.533, -0.3127, -0.501, -0.1179, -0.5257, -0.1527,
        0.208, 0.6069, 0.4717, 0.6916, 0.3173, 0.3884, 0.1149, -0.2679,
        -0.2679, -0.2679
    ]
    assert b_list == approx3(b_list_ref)


def test_aggrescan3d():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    psico.aggrescanning.aggrescan3d()
    b_list = psico.querying.iterate_to_list("guide", "b")
    b_list_ref = [
        -3.0073, -3.3591, 0.0, -0.9979, 0.0, -1.757, -0.3646, 1.2682, 0.198,
        0.0, -1.9834, -1.2031, 0.0, -1.1033, 0.0, -4.0497, 0.0, -4.1411,
        -2.668, -2.9765, -3.1231, -4.0287, 0.0, -4.5619, -2.1746, 0.0, 0.0,
        -1.7524, -3.1351, 0.0, -2.6589, -4.1382, -4.4392, -3.0258, 0.0,
        -0.6704, 0.2097, -1.268, -2.9283, -2.3215, 0.0, -3.0799, 0.0, -1.1319,
        0.1755, 0.2896, 0.0, -1.9384, -2.6352, 0.0, -5.6244, -5.226, 0.0,
        -6.2019, -4.6744, 0.0, -3.675, -4.1248, -2.6781, -2.1701, 0.0, -3.4446,
        -4.9271, -4.4772, -2.2652, -1.0839, 0.0, -0.6793, 0.0, -0.2731,
        -0.8327, -3.656, -1.8123, -4.2951, 0.0, 0.0
    ]
    assert b_list == approx3(b_list_ref)
