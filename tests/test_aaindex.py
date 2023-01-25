import psico.aaindex
import psico.querying
from pymol import cmd
from pathlib import Path
from pytest import approx

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_hydropathy2b(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')

    # use cached sample file, skip download
    cmd.set("fetch_path", str(DATA_PATH / "aaindex"))

    psico.aaindex.hydropathy2b()

    b_list = psico.querying.iterate_to_list("guide", "b")
    assert b_list == approx([3.8, -0.8, 1.9, 1.9, -1.6, -3.5])


def test_hydropathy2b_exposed(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')
    cmd.flag("ignore", "resn MSE", "clear")

    # use cached sample file, skip download
    cmd.set("fetch_path", str(DATA_PATH / "aaindex"))

    psico.aaindex.hydropathy2b_exposed()

    b_list = psico.querying.iterate_to_list("guide", "b")
    assert b_list == approx([3.8, -0.8, 1.9, 1.9, -1.6, -3.5])

    q_list = psico.querying.iterate_to_list("guide", "q")
    assert q_list == approx([0.76, 0.97, 0.86, 0.79, 0.88, 0.89], abs=1e-2)

    color_list = psico.querying.iterate_to_list("guide", "color")
    color_list = [cmd.get_color_tuple(color) for color in color_list]
    assert color_list == [
        approx([0.2, 0.6, 0.2]),
        approx([1.0, 1.0, 1.0]),
        approx([0.5412, 0.7686, 0.5412], abs=1e-3),
        approx([0.5804, 0.7882, 0.5804], abs=1e-3),
        approx([1.0, 1.0, 1.0]),
        approx([1.0, 1.0, 1.0]),
    ]

    psico.aaindex.hydropathy2b_exposed(selection="all",
                                       state=1,
                                       palette="0xFF0000 0x0000FF")

    color_list = psico.querying.iterate_to_list("guide", "color")
    color_list = [cmd.get_color_tuple(color) for color in color_list]
    assert color_list == [
        approx([0.0, 0.0, 1.0]),
        approx([1.0, 0.0, 0.0]),
        approx([0.4274, 0.0, 0.5686], abs=1e-3),
        approx([0.4745, 0.0, 0.5215], abs=1e-3),
        approx([1.0, 0.0, 0.0]),
        approx([1.0, 0.0, 0.0]),
    ]


def test_aaindex2b(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')

    # use cached sample file, skip download
    cmd.set("fetch_path", str(DATA_PATH / "aaindex"))

    psico.aaindex.aaindex2b("LAWE840101", "name CA", "q")

    b_list = psico.querying.iterate_to_list("guide", "q")
    assert b_list == approx([1.02, 0.05, 0.81, 0.81, 2.03, -0.75])
