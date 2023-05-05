import psico.setting
from pymol import cmd


def test_save_settings(tmp_path):
    cmd.reinitialize()
    filepath = tmp_path / "mysettings.py"
    cmd.set("fetch_path", "abc")
    psico.setting.save_settings(str(filepath))
    content = filepath.read_text()
    assert "fetch_path" in content
    assert "abc" in content


def test_paper_settings():
    cmd.reinitialize()
    assert cmd.get_setting_int("cartoon_side_chain_helper") < 1
    psico.setting.paper_settings()
    assert cmd.get_setting_int("cartoon_side_chain_helper") == 1
    assert cmd.get_setting_int("cartoon_fancy_helices") < 1
    assert cmd.get_setting_int("ray_shadows") == 0
    psico.setting.paper_settings(fancy=1)
    assert cmd.get_setting_int("cartoon_fancy_helices") == 1


def test_set_temporary():
    cmd.reinitialize()
    namea = "cartoon_sampling"
    nameb = "cartoon_gap_cutoff"
    cmd.set(namea, 1)
    cmd.set(nameb, 5)
    assert cmd.get_setting_int(namea) == 1
    assert cmd.get_setting_int(nameb) == 5
    with psico.setting.set_temporary((namea, 2), (nameb, 7)):
        assert cmd.get_setting_int(namea) == 2
        assert cmd.get_setting_int(nameb) == 7
        with psico.setting.set_temporary(**{namea: 3, nameb: 9}):
            assert cmd.get_setting_int(namea) == 3
            assert cmd.get_setting_int(nameb) == 9
        assert cmd.get_setting_int(namea) == 2
        assert cmd.get_setting_int(nameb) == 7
    assert cmd.get_setting_int(namea) == 1
    assert cmd.get_setting_int(nameb) == 5
