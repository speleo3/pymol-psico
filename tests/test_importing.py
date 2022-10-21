import psico.importing
from pytest import approx
from pymol import cmd
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / 'data'

ALA_TRAJ_COORDS_STATE_3 = [
    [-0.669, -1.228, -0.491],
    [0.004, 0.059, -0.495],
    [-0.502, 0.844, 0.715],
    [1.493, -0.116, -0.492],
    [2.018, -1.228, -0.509],
    [-1.483, 0.49, 1.096],
    [0.173, 0.82, 1.6],
    [-0.651, 1.917, 0.491],
    [-0.123, -2.152, -0.496],
    [-0.266, 0.597, -1.411],
]

ALA_TRAJ_COORDS_STATE_8 = [
    [-0.667, -1.194, -0.453],
    [0.002, 0.091, -0.494],
    [-0.508, 0.86, 0.729],
    [1.489, -0.072, -0.536],
    [2.278, -0.728, 0.133],
    [-1.004, 0.202, 1.462],
    [0.312, 1.376, 1.256],
    [-1.245, 1.627, 0.438],
    [-1.042, -2.13, -0.82],
    [-0.268, 0.608, -1.424],
]


def test_loadall_traj():
    cmd.reinitialize()
    cmd.set("retain_order")
    cmd.fragment("ala", "m1")
    assert cmd.count_states("m1") == 1
    psico.importing.loadall_traj(DATA_PATH / "ala-traj-*.dcd")
    assert cmd.count_states("m1") == 10
    psico.importing.loadall_traj(DATA_PATH / "ala-traj-*.dcd", state="append")
    assert cmd.count_states("m1") == 20
    i = 8
    coords = lambda state: cmd.get_coords("m1", state=state).tolist()[i]
    assert coords(3) == approx(ALA_TRAJ_COORDS_STATE_3[i], abs=1e-3)
    assert coords(8) == approx(ALA_TRAJ_COORDS_STATE_8[i], abs=1e-3)
    assert coords(13) == approx(ALA_TRAJ_COORDS_STATE_3[i], abs=1e-3)
    assert coords(18) == approx(ALA_TRAJ_COORDS_STATE_8[i], abs=1e-3)
