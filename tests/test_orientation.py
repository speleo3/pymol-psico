import psico.orientation
from pathlib import Path
from pymol import cmd
from pytest import approx

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_angle_between_helices():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1nmr-frag-nohydro.pdb", "m1")

    angle = psico.orientation.angle_between_helices("resi 53-61",
                                                    "resi 63-79",
                                                    method="helix")
    assert angle == approx(159.164, abs=0.01)
    assert set(["oriVec01", "oriVec02", "angle01"]).issubset(cmd.get_names())

    angle = psico.orientation.angle_between_helices("resi 53-61",
                                                    "resi 63-79",
                                                    method="cafit")
    assert angle == approx(158.624, abs=0.01)
    assert set(["oriVec03", "oriVec04", "angle02"]).issubset(cmd.get_names())


# def loop_orientation(selection, state=STATE, visualize=1, quiet=1):
# def plane_orientation(selection, state=STATE, visualize=1, guide=0, quiet=1):
# def angle_between_domains(selection1, selection2, method='align',
