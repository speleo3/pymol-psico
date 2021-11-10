'''
(c) 2011 Thomas Holder, MPI for Develpomental Biology

License: BSD-2-Clause
'''

from pymol import cmd

def qdelaunay(X, n=0, m=0, options='Qt', qdelaunay_exe='qdelaunay'):
    '''
Triangulation using qdelaunay. (http://www.qhull.org)

@param X: iterable object of points (not numpy.matrix). If n or m are 0, then
          X must be sequence.
@returns: iterator of regions
    '''
    import subprocess
    if not n: n = len(X[0])
    if not m: m = len(X)
    process = subprocess.Popen([qdelaunay_exe, 'i'] + options.split(),
            universal_newlines=True,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    # input
    print(n, file=process.stdin)
    print(m, file=process.stdin)
    for coord in X:
        print('%f %f %f' % tuple(coord), file=process.stdin)
    process.stdin.close()
    # output
    out_it = iter(process.stdout)
    n_regions = next(out_it)
    for line in out_it:
        a = line.split()
        assert len(a) == n + 1, 'Wrong number of vertices in output'
        yield list(map(int, a))

def delaunay(selection='enabled', name=None, cutoff=10.0, as_cgo=0,
        qdelaunay_exe='qdelaunay', state=-1, quiet=1):
    '''
DESCRIPTION

    Full-atom Delaunay Tessalator

    Creates either a molecular object with delaunay edges as bonds, or a CGO
    object with edge colors according to edge length.

USAGE

    delaunay [ selection [, name [, cutoff=10.0 [, as_cgo=0 ]]]]

SEE ALSO

    PyDeT plugin: http://pymolwiki.org/index.php/PyDet
    '''
    from chempy import cpv, Bond
    if name is None:
        name = cmd.get_unused_name('delaunay')
    cutoff = float(cutoff)
    as_cgo = int(as_cgo)
    state, quiet = int(state), int(quiet)
    if state < 1:
        state = cmd.get_state()
    model = cmd.get_model(selection, state)
    regions_iter = qdelaunay((a.coord for a in model.atom), 3, len(model.atom),
            qdelaunay_exe=qdelaunay_exe)
    edges = set(tuple(sorted([region[i-1], region[i]]))
            for region in regions_iter for i in range(len(region)))

    edgelist=[]
    r = []

    minco = 9999
    maxco = 0

    for edge in edges:
        ii, jj = edge
        a = model.atom[ii]
        b = model.atom[jj]
        co = cpv.distance(a.coord, b.coord)
        if cutoff > 0.0 and co > cutoff:
            continue
        if as_cgo:
            minco=min(co,minco)
            maxco=max(co,maxco)
            edgelist.append(a.coord + b.coord + [co])
        else:
            bnd = Bond()
            bnd.index = [ii, jj]
            model.add_bond(bnd)
        r.append((a,b,co))

    if not as_cgo:
        cmd.load_model(model, name, 1)
        return r

    from pymol.cgo import CYLINDER

    difco = maxco-minco
    obj = []
    mm = lambda x: max(min(x, 1.0), 0.0)
    for e in edgelist:
        co = ((e[6]-minco)/difco)**(0.75)
        color = [mm(1-2*co), mm(1-abs(2*co-1)), mm(2*co-1)]
        obj.extend([CYLINDER] + e[0:6] + [0.05] + color + color)

    cmd.load_cgo(obj, name)
    return r

cmd.extend('delaunay', delaunay)

# tab-completion of arguments
cmd.auto_arg[0]['delaunay'] = cmd.auto_arg[0]['zoom']

# vi:sw=4:expandtab:smarttab
