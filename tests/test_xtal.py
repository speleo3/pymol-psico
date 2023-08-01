import psico.xtal
from pathlib import Path
from pymol import cmd
from pytest import approx

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_cellbasis():
    cmd.reinitialize()
    cellbasis = psico.xtal.cellbasis([90, 90, 90], [12.3, 45.6, 78.9])
    assert cellbasis.tolist() == [
        approx([12.3, 0, 0, 0]),
        approx([0., 45.6, 0., 0.]),
        approx([0., 0., 78.9, 0.]),
        approx([0., 0., 0., 1.]),
    ]


COORDS_MATES_1UBQ_20 = [
    approx(dim, abs=1e-3) for dim in [
        [
            27.361, 28.054, 29.258, 29.93, 28.523, 28.946, 48.899, 48.206,
            47.002, 46.33, 47.737, 47.314, 1.941, 2.634, 3.838, 4.51, 3.103,
            3.526, 23.479, 22.786, 21.582, 20.91, 22.317, 21.894
        ],
        [
            17.959, 16.835, 17.318, 16.477, 15.82, 16.445, 24.811, 25.935,
            25.452, 26.293, 26.95, 26.325, 46.196, 47.32, 46.837, 47.678,
            48.335, 47.71, -3.426, -4.55, -4.067, -4.908, -5.565, -4.94
        ],
        [
            8.559, 9.21, 9.984, 10.606, 8.182, 6.967, -5.916, -5.265, -4.491,
            -3.869, -6.293, -7.508, 20.391, 19.74, 18.966, 18.344, 20.768,
            21.983, 34.866, 34.215, 33.441, 32.819, 35.243, 36.458
        ],
    ]
]


def test_supercell():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    psico.xtal.supercell()
    assert cmd.get_names() == [
        "m1",
        "m000_1",
        "m000_2",
        "m000_3",
        "m000_4",
        "m000",
        "supercell",
    ]
    assert cmd.count_atoms("m000_1") == 660
    assert cmd.count_atoms("m000_4") == 660
    assert cmd.get_type("m000") == "object:group"
    assert cmd.get_type("supercell") == "object:cgo"
    coords = cmd.get_coords("m000_* & resi 20").round(5).T.tolist()
    assert coords == COORDS_MATES_1UBQ_20
    assert cmd.get_object_matrix("m000_1") == approx([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    assert cmd.get_object_matrix("m000_2") == approx([-1, 0, 0, 76.26, 0, -1, 0, 42.77, 0, 0, 1, -14.475, 0, 0, 0, 1])
    assert cmd.get_object_matrix("m000_3") == approx([1, 0, 0, -25.42, 0, -1, 0, 64.155, 0, 0, -1, 28.95, 0, 0, 0, 1])
    assert cmd.get_object_matrix("m000_4") == approx([-1, 0, 0, 50.84, 0, 1, 0, -21.385, 0, 0, -1, 43.425, 0, 0, 0, 1])


def test_pdbremarks(tmp_path):
    path = tmp_path / "foo.pdb"
    path.write_text("""HEADER    Test
REMARK   1 """ """
REMARK   1 REFERENCE 1
REMARK   1  AUTH   M.MUSTERMANN
REMARK   1  TITL   TESTING CODE
REMARK 100     """ """
REMARK 100 THIS ENTRY HAS BEEN PROCESSED MANUALLY.
""")
    remarks = psico.xtal.pdbremarks(str(path))
    assert remarks == {
        1: ["\n", "REFERENCE 1\n", " AUTH   M.MUSTERMANN\n", " TITL   TESTING CODE\n"],
        100: ["    \n", "THIS ENTRY HAS BEEN PROCESSED MANUALLY.\n"],
    }


# def biomolecule(name=None, filename=None, prefix=None, number=1, suffix=None,


def test_cell_axes():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    psico.xtal.cell_axes(name="a1")
    assert cmd.get_type("a1") == "object:cgo"
