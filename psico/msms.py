'''
(c) 2016 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def save_xyzr(filename, selection='all', state=1, _colorsout=None, *, _self=cmd):
    '''
DESCRIPTION

    Write the given selection to an xyzr or xyzrn (determined by extension)
    file for MSMS.
    '''
    if filename.endswith('xyzrn'):
        expr = 'callback(x, y, z, vdw, name, resn, resi)'
        fmt = '%.3f %.3f %.3f %.2f 1 %s_%s_%s\n'
    else:
        expr = 'callback(x, y, z, vdw)'
        fmt = '%.3f %.3f %.3f %.2f\n'

    if isinstance(_colorsout, list):
        expr += ';_colorsout.append(color)'

    handle = open(filename, 'w')
    try:
        _self.iterate_state(state, selection, expr, space={
            'callback': lambda *args: handle.write(fmt % args),
            '_colorsout': _colorsout})
    finally:
        handle.close()


def load_msms_surface(filename, name='', _colors=None, *, _self=cmd):
    '''
DESCRIPTION

    Load MSMS .vert and .face files as a CGO
    '''
    from pymol import cgo
    from pymol.cgo import NORMAL, VERTEX, COLOR

    if _colors:
        _colors = [_self.get_color_tuple(c) for c in _colors]

    if filename.endswith('.vert') or filename.endswith('.face'):
        filename = filename[:-5]

    # vertex file
    line_iter = iter(open(filename + '.vert'))

    # skip header
    for line in line_iter:
        if not line.startswith('#'):
            break

    # read vertices
    vertices = [None]  # make 1-indexable
    for line in line_iter:
        data = line.split()
        vertex = [float(x) for x in data[0:3]]
        normal = [float(x) for x in data[3:6]]
        sphere = int(data[7]) - 1
        vertices.append((vertex, normal, sphere))

    # faces file
    line_iter = iter(open(filename + '.face'))

    # skip header
    for line in line_iter:
        if not line.startswith('#'):
            break

    cgobuf = [cgo.BEGIN, cgo.TRIANGLES]

    # read triangles
    for line in line_iter:
        for index in line.split()[:3]:
            data = vertices[int(index)]
            if _colors:
                cgobuf.append(COLOR)
                cgobuf.extend(_colors[data[2]])
            cgobuf.append(NORMAL)
            cgobuf.extend(data[1])
            cgobuf.append(VERTEX)
            cgobuf.extend(data[0])

    cgobuf.append(cgo.END)

    if not name:
        name = _self.get_unused_name('msmssurf')

    _self.load_cgo(cgobuf, name)


def _get_solvent_radius(selection, state, *, _self=cmd):
    '''Get object-state level solvent_radius'''
    radius = [None]

    try:
        _self.iterate_state(state,
            'first ({})'.format(selection),
            'radius[0] = s.solvent_radius',
            space=locals())

        if radius[0] is not None:
            # Note: One of my test cases failed with radius 2.75, rounding
            # to one digit worked
            return '{:.1}'.format(radius[0])
    except:
        print('Using global solvent_radius')

    return _self.get('solvent_radius')


def msms_surface(selection='polymer', state=1, density=3, name='',
        atomcolors=0, exe='msms', preserve=0, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Run MSMS on the given selection and load the generated surface
    as a CGO.

    Note that PyMOL's own surface generation by default excludes atoms with
    the "ignore" flag. To match this behavior, use "polymer" (recommended)
    or "not flag 25" in your selection.

ARGUMENTS

    selection = str: atom selection {default: polymer}

    state = int: object state {default: 1}

    density = float: MSMS surface point density {default: 3}

    name = str: name of CGO object to create

    atomcolors = 0/1: color surface by atom colors {default: 0}

EXAMPLE

    # optional: use Connolly radii
    atmtypenumbers /tmp/MSMS-release/atmtypenumbers

    # show surface for protein
    msms_surface polymer
    '''
    import os, tempfile, subprocess, shutil

    density = float(density)
    hdensity = density + 2

    colors = [] if int(atomcolors) else None

    tmpdir = tempfile.mkdtemp()
    tmp_if = os.path.join(tmpdir, 'xxxx.xyzr')
    tmp_of = os.path.join(tmpdir, 'xxxx')

    try:
        save_xyzr(tmp_if, selection, state, _colorsout=colors, _self=_self)

        subprocess.check_call([exe,
            '-density', str(density),
            '-hdensity', str(hdensity),
            '-if', tmp_if,
            '-of', tmp_of,
            '-no_area',
            '-probe_radius', _get_solvent_radius(selection, state, _self=_self),
        ], cwd=tmpdir)

        load_msms_surface(tmp_of, name, colors, _self=_self)

    except OSError:
        raise CmdException('Cannot execute exe=' + exe)
    finally:
        if not int(preserve):
            shutil.rmtree(tmpdir)
        elif not int(quiet):
            print(' Notice: not deleting ' + tmpdir)


def atmtypenumbers(filename='atmtypenumbers', selection='all', united=1,
        quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Update VDW radii for the given selection, based on "atmtypenumbers" file.

ARGUMENTS

    filename = str: path to "atmtypenumbers" file

    selection = str: atom selection to update {default: all}

    united = 1: Use united (implicit hydrogens) radii {default}
    united = 0: Use explicit radii

EXAMPLE

    atmtypenumbers /tmp/MSMS-release/atmtypenumbers, united=0
    '''
    import re

    united = int(united)
    quiet = int(quiet)

    if (united and not quiet and
            _self.count_atoms('(%s) and hydro' % (selection))):
        print(" Warning: united=1 but found hydrogens in selection")

    types = {}
    patterns = []

    def _true_func(s):
        return True

    def _get_match_func(p):
        if p == '*':
            return _true_func
        return re.compile(p.replace('_', ' ') + '$').match

    try:
        handle = open(filename)
    except IOError:
        raise CmdException("can't open file '%s', please provide correct "
                "filename" % (filename))

    for line in handle:
        fields = line.split()
        for i, field in enumerate(fields):
            if field.startswith('#'):
                fields = fields[:i]
                break
        if not fields:
            continue

        if fields[0] == 'radius':
            vdwidx = 4 if united and len(fields) > 4 else 3
            types[fields[1]] = float(fields[vdwidx])
            continue

        patterns.append((
            _get_match_func(fields[0]),
            _get_match_func(fields[1]),
            types[fields[2]]))

    handle.close()

    def callback(resn, name, vdw):
        for p in patterns:
            if p[0](resn) and p[1](name):
                return p[2]
        if not quiet:
            print(" Warning: no match for '%s/%s'" % (resn, name))
        return vdw

    _self.alter(selection, 'vdw = callback(resn, name, vdw)',
            space={'callback': callback})
    _self.rebuild(selection)


cmd.extend('save_xyzr', save_xyzr)
cmd.extend('load_msms_surface', load_msms_surface)
cmd.extend('msms_surface', msms_surface)
cmd.extend('atmtypenumbers', atmtypenumbers)

# auto-completion
cmd.auto_arg[0]['msms_surface'] = cmd.auto_arg[1]['select']
