import psico.numeric
import pytest


def test_pdist_squareform():
    try:
        import scipy.spatial.distance  # noqa: F401
    except ModuleNotFoundError:
        pytest.skip("scipy not available")

    import numpy.random
    X = numpy.random.random((5, 3)) * 100

    for method in ("euclidean", "sqeuclidean"):
        D = psico.numeric.pdist_squareform(X, method)
        D_numpy = psico.numeric.pdist_squareform_numpy(X, method)
        D_scipy = psico.numeric.pdist_squareform_scipy(X, method)
        assert numpy.allclose(D, D_numpy, atol=1e-4)
        assert numpy.allclose(D, D_scipy, atol=1e-4)


def test_pdist_squareform_numpy():
    import numpy
    X = numpy.array([
        [0, 0, 0],
        [3, 0, 0],
    ])
    D = psico.numeric.pdist_squareform_numpy(X)
    assert numpy.allclose(D, [
        [0, 3],
        [3, 0],
    ])
    D = psico.numeric.pdist_squareform_numpy(X, "sqeuclidean")
    assert numpy.allclose(D, [
        [0, 9],
        [9, 0],
    ])
