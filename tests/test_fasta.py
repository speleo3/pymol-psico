import psico.fasta
from pymol import cmd


def test_fasta(capsys, tmp_path):
    cmd.reinitialize()
    cmd.fab("A" * 20 + "CDEF", "m1", chain="B")
    cmd.fab("X/1/ EFG Y/1/ HIK", "m2")
    cmd.alter("m1 & resi 11-", "resv += 5")
    psico.fasta.fasta(wrap=20, gapped=0)
    captured = capsys.readouterr()
    assert captured.out.split() == [
        ">m1_B", "A" * 20, "CDEF", ">m2_X", "EFG", ">m2_Y", "HIK"
    ]
    psico.fasta.fasta("m1")
    captured = capsys.readouterr()
    assert captured.out.split() == [
        ">m1_B", "A" * 10 + "-" * 5 + "A" * 10 + "CDEF"
    ]
    filepath = tmp_path / "out.fasta"
    psico.fasta.fasta(filename=filepath)
    assert filepath.is_file()


def test_pir(capsys):
    cmd.reinitialize()
    cmd.fab("ACD", "m1")
    psico.fasta.pir()
    captured = capsys.readouterr()
    assert captured.out.split() == [">P1;m1", "structure:m1:1::3:::::", "ACD*"]


def test_save_colored_fasta(tmp_path):
    cmd.reinitialize()
    cmd.fab("ACD", "m1")
    filepath = tmp_path / "out.fasta"
    psico.fasta.save_colored_fasta(filepath)
    assert filepath.is_file()
