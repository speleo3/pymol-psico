import psico.moving
from pymol import cmd

# def frames2states(selection, specification):
# def save_movie_mpeg1(filename, mode='', first=0, last=0, preserve=0,
# def matrix_to_ttt(names, reverse=0, state=-1, quiet=1):


def test_get_keyframes():
    cmd.reinitialize()
    assert psico.moving.get_keyframes() == []
    cmd.mset("1x30")
    cmd.mview("store", 4)
    cmd.mview("store", 22)
    assert psico.moving.get_keyframes() == [4, 22]


def test_closest_keyframe():
    cmd.reinitialize()
    assert psico.moving.closest_keyframe() is None
    cmd.mset("1x30")
    cmd.mview("store", 4)
    cmd.mview("store", 22)
    assert psico.moving.closest_keyframe() == 4
    assert cmd.get_frame() == 4
    cmd.frame(20)
    assert psico.moving.closest_keyframe() == 22
    assert cmd.get_frame() == 22


def test_next_keyframe():
    cmd.reinitialize()
    cmd.mset("1x30")
    cmd.mview("store", 4)
    cmd.mview("store", 9)
    cmd.mview("store", 22)
    assert psico.moving.next_keyframe() == 4
    assert psico.moving.next_keyframe() == 9


def test_prev_keyframe():
    cmd.reinitialize()
    cmd.mset("1x30")
    cmd.mview("store", 4)
    cmd.mview("store", 9)
    cmd.mview("store", 22)
    assert psico.moving.prev_keyframe() is None
    cmd.frame(12)
    assert psico.moving.prev_keyframe() == 9


def test_get_mdo_commands():
    cmd.reinitialize()
    cmd.mset("1x6")
    cmd.mdo(3, "color blue")
    cmd.mappend(3, "bg yellow")
    cmd.mdo(5, "show sticks")
    assert psico.moving.get_mdo_commands() == [
        "", "", "color blue;bg yellow", "", "show sticks", ""
    ]


# dump_mviews
