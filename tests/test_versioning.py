from psico.versioning import make_version_int_tuple


def test_make_version_int_tuple():
    assert make_version_int_tuple("0") == ()
    assert make_version_int_tuple("0.99rc6") == (0, 99)
    assert make_version_int_tuple("1") == (1,)
    assert make_version_int_tuple("1.0") == (1,)
    assert make_version_int_tuple("1.1") == (1, 1)
    assert make_version_int_tuple("1.2r2") == (1, 2)
    assert make_version_int_tuple("2.6.0a0") == (2, 6)
    assert make_version_int_tuple("2.6.0") == (2, 6)
    assert make_version_int_tuple("2.6.3") == (2, 6, 3)
    assert make_version_int_tuple(".") == ()
    assert make_version_int_tuple("foo") == ()
    assert make_version_int_tuple("foo.300") == ()
    assert make_version_int_tuple("300") == (300,)
    assert make_version_int_tuple("300") == (300,)
