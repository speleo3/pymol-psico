'''
numpy and scipy helpers.

(c) 2021 Thomas Holder

License: BSD-2-Clause
'''


def pdist_squareform_numpy(X, metric="euclidean"):
    """Distance matrix. Uses numpy.
    """
    import numpy

    X = numpy.asarray(X)
    Q = (X**2).sum(1)
    S = numpy.add.outer(Q, Q)
    D = (S - 2 * (X @ X.T)).clip(0)  # sqeuclidean

    if metric == "sqeuclidean":
        return D
    elif metric == "euclidean":
        return D**0.5

    raise ValueError(f"Unknown Distance Metric: {metric}")


def pdist_squareform_scipy(X, metric="euclidean"):
    """Distance matrix. Uses scipy.
    """
    from scipy.spatial.distance import pdist, squareform
    return squareform(pdist(X, metric))


def pdist_squareform(X, metric="euclidean"):
    """Distance matrix. Uses scipy if available (faster) or numpy.
    """
    try:
        return pdist_squareform_scipy(X, metric)
    except ModuleNotFoundError:
        return pdist_squareform_numpy(X, metric)
