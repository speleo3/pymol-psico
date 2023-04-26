import psico.helping
from pymol import cmd


def test_grepset(capsys):
    cmd.reinitialize()
    psico.helping.grepset('^ray_shadow$')
    captured = capsys.readouterr()
    assert captured.out.startswith("ray_shadow : on")
    assert "1 settings matched" in captured.out
    psico.helping.grepset('ray_shadow')
    captured = capsys.readouterr()
    assert "4 settings matched" in captured.out


def test_apropos(capsys):
    cmd.reinitialize()
    psico.helping.apropos('load')
    captured = capsys.readouterr()
    assert "EXACT MATCH FOR: load" in captured.out
    assert captured.out.count(" : ") > 10
    psico.helping.apropos('ABC some unfindable string XYZ')
    captured = capsys.readouterr()
    assert captured.out == ""


def test_api_info(capsys):
    cmd.reinitialize()
    func = psico.helping.api_info("ray")
    captured = capsys.readouterr()
    assert " CMD: ray" in captured.out
    assert " API: pymol.viewing.ray" in captured.out
    assert " FILE: " in captured.out
    assert func is cmd.ray


def test_write_html_ref(tmp_path):
    path = tmp_path / "ref.html"
    psico.helping.write_html_ref(path, prefix="pymol.viewing")
    content = path.read_text()
    assert content.count("<hr ") > 5


def test_write_txt_ref(tmp_path):
    path = tmp_path / "ref.txt"
    psico.helping.write_txt_ref(path, prefix="pymol.viewing")
    content = path.read_text()
    assert '"zoom" scales and translates' in content
