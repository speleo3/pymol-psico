'''
Simplified morphing workflow

(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd


def morpheasy(source, target, source_state=0, target_state=0, name=None,
        refinement=5, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Morph source to target, based on sequence alignment

USAGE

    morpheasy source, target [, source_state [, target_state [, name ]]]

EXAMPLE

    fetch 1akeA 4akeA, async=0
    extra_fit
    morpheasy 1akeA, 4akeA
    '''
    try:
        from epymol import rigimol
    except ImportError:
        print('No epymol available, please use a "Incentive PyMOL" build')
        print('You may use "morpheasy_linear" instead')
        return

    from .editing import mse2met
    from .querying import get_selection_state

    # arguments
    source_state = int(source_state)
    target_state = int(target_state)
    refinement = int(refinement)
    quiet = int(quiet)

    if source_state < 1:
        source_state = get_selection_state(source, _self=_self)
    if target_state < 1:
        target_state = get_selection_state(target, _self=_self)

    # temporary objects
    # IMPORTANT: cmd.get_raw_alignment does not work with underscore object names!
    alnobj = _self.get_unused_name('_aln')
    so_obj = _self.get_unused_name('source')  # see above
    ta_obj = _self.get_unused_name('target')  # see above
    so_sel = _self.get_unused_name('_source_sel')
    ta_sel = _self.get_unused_name('_target_sel')
    _self.create(so_obj, source, source_state, 1)
    _self.create(ta_obj, target, target_state, 1)
    mse2met(so_obj, _self=_self)
    mse2met(ta_obj, _self=_self)

    # align sequence
    _self.align(ta_obj, so_obj, object=alnobj, cycles=0, transform=0,
            mobile_state=1, target_state=1)
    _self.refresh()
    _self.select(so_sel, '%s and %s' % (so_obj, alnobj))
    _self.select(ta_sel, '%s and %s' % (ta_obj, alnobj))
    alnmap = dict(_self.get_raw_alignment(alnobj))
    alnmap.update(dict((v, k) for (k, v) in alnmap.items()))

    # copy source atom identifiers to temporary target
    idmap = dict()
    _self.iterate(so_sel, 'idmap[model,index] = (segi,chain,resi,resn,name)',
            space={'idmap': idmap})
    _self.alter(ta_sel, '(segi,chain,resi,resn,name) = idmap[alnmap[model,index]]',
            space={'idmap': idmap, 'alnmap': alnmap})

    # remove unaligned
    _self.remove('%s and not %s' % (so_obj, so_sel))
    _self.remove('%s and not %s' % (ta_obj, ta_sel))
    assert _self.count_atoms(so_obj) == _self.count_atoms(ta_obj)
    _self.sort(so_obj)
    _self.sort(ta_obj)

    # append target to source as 2-state morph-in object
    _self.create(so_obj, ta_obj, 1, 2)

    # morph
    if name is None:
        name = _self.get_unused_name('morph')
    rigimol.morph(so_obj, name, refinement=refinement, **{'async': 0})

    # clean up
    for obj in [alnobj, so_obj, so_sel, ta_obj, ta_sel]:
        _self.delete(obj)

    return name


def morpheasy_linear(source, target, source_state=0, target_state=0, name=None,
        steps=30, match='align', quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Morph by linear interpolation in cartesian space (like LSQMAN).

    This is the poor man's version of morphing, it's quick but will produce
    distorted intermediate conformations. Does not require rigimol (incentive
    PyMOL product). Requires numpy.

SEE ALSO

    morpheasy
    '''
    from numpy import asfarray
    from .fitting import MatchMaker
    from .querying import get_selection_state

    # arguments
    source_state = int(source_state)
    target_state = int(target_state)
    steps, quiet = int(steps), int(quiet)

    if source_state < 1:
        source_state = get_selection_state(source, _self=_self)
    if target_state < 1:
        target_state = get_selection_state(target, _self=_self)

    mm = MatchMaker(source, target, match, _self=_self)

    csource = asfarray(_self.get_coords(mm.mobile, source_state))
    ctarget = asfarray(_self.get_coords(mm.target, target_state))
    cdiff = ctarget - csource

    if name is None:
        name = _self.get_unused_name('morph')
    _self.create(name, mm.mobile, source_state, 1)

    for state in range(2, steps + 1):
        c = csource + cdiff * float(state - 1) / (steps - 1)
        _self.load_coordset(c.tolist(), name, state)

    return name


cmd.extend('morpheasy', morpheasy)
cmd.extend('morpheasy_linear', morpheasy_linear)

cmd.auto_arg[0].update([
    ('morpheasy', cmd.auto_arg[0]['align']),
    ('morpheasy_linear', cmd.auto_arg[0]['align']),
])
cmd.auto_arg[1].update([
    ('morpheasy', cmd.auto_arg[1]['align']),
    ('morpheasy_linear', cmd.auto_arg[1]['align']),
])

# vi: ts=4:sw=4:smarttab:expandtab
