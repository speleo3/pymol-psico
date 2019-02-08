'''
Simplified morphing workflow

(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd

def morpheasy(source, target, source_state=0, target_state=0, name=None,
        refinement=5, quiet=1):
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

    if source_state < 1: source_state = get_selection_state(source)
    if target_state < 1: target_state = get_selection_state(target)

    # temporary objects
    # IMPORTANT: cmd.get_raw_alignment does not work with underscore object names!
    alnobj = cmd.get_unused_name('_aln')
    so_obj = cmd.get_unused_name('source') # see above
    ta_obj = cmd.get_unused_name('target') # see above
    so_sel = cmd.get_unused_name('_source_sel')
    ta_sel = cmd.get_unused_name('_target_sel')
    cmd.create(so_obj, source, source_state, 1)
    cmd.create(ta_obj, target, target_state, 1)
    mse2met(so_obj)
    mse2met(ta_obj)

    # align sequence
    cmd.align(ta_obj, so_obj, object=alnobj, cycles=0, transform=0,
            mobile_state=1, target_state=1)
    cmd.refresh()
    cmd.select(so_sel, '%s and %s' % (so_obj, alnobj))
    cmd.select(ta_sel, '%s and %s' % (ta_obj, alnobj))
    alnmap = dict(cmd.get_raw_alignment(alnobj))
    alnmap.update(dict((v,k) for (k,v) in alnmap.items()))

    # copy source atom identifiers to temporary target
    idmap = dict()
    cmd.iterate(so_sel, 'idmap[model,index] = (segi,chain,resi,resn,name)',
            space={'idmap': idmap})
    cmd.alter(ta_sel, '(segi,chain,resi,resn,name) = idmap[alnmap[model,index]]',
            space={'idmap': idmap, 'alnmap': alnmap})

    # remove unaligned
    cmd.remove('%s and not %s' % (so_obj, so_sel))
    cmd.remove('%s and not %s' % (ta_obj, ta_sel))
    assert cmd.count_atoms(so_obj) == cmd.count_atoms(ta_obj)
    cmd.sort(so_obj)
    cmd.sort(ta_obj)

    # append target to source as 2-state morph-in object
    cmd.create(so_obj, ta_obj, 1, 2)

    # morph
    if name is None:
        name = cmd.get_unused_name('morph')
    rigimol.morph(so_obj, name, refinement=refinement, **{'async': 0})

    # clean up
    for obj in [alnobj, so_obj, so_sel, ta_obj, ta_sel]:
        cmd.delete(obj)

    return name

def morpheasy_linear(source, target, source_state=0, target_state=0, name=None,
        steps=30, match='align', quiet=1):
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
    from .fitting import matchmaker
    from .importing import load_coords
    from .querying import get_selection_state, get_coords

    # arguments
    source_state = int(source_state)
    target_state = int(target_state)
    steps, quiet = int(steps), int(quiet)

    if source_state < 1: source_state = get_selection_state(source)
    if target_state < 1: target_state = get_selection_state(target)

    msource, mtarget, tmp_names = matchmaker(source, target, match)

    csource = asfarray(get_coords(msource, source_state))
    ctarget = asfarray(get_coords(mtarget, target_state))
    cdiff = ctarget - csource

    if name is None:
        name = cmd.get_unused_name('morph')
    cmd.create(name, msource, source_state, 1)

    for state in range(2, steps+1):
        c = csource + cdiff * float(state-1) / (steps-1)
        load_coords(c.tolist(), name, state)

    # clean up
    for obj in tmp_names:
        cmd.delete(obj)

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
