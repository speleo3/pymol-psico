'''
Viewing and Coloring Stuff

(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def nice(selection='(all)', simple=1e5, *, _self=cmd):
    '''
DESCRIPTION

    My favourite representation

    Color: by chain (elem C) and by element (others)

    Representation: cartoon (polymer), sticks (organic),
    spheres (inorganic) and nonbonded (solvent)

USAGE

    nice [ selection ]
    '''
    simple = int(simple)
    _self.util.cbc(selection)
    _self.color('atomic', '(%s) and not elem C' % (selection))
    if simple and _self.count_atoms(selection) >= simple:
        _self.show_as('ribbon', selection)
        _self.show_as('lines', '(%s) and organic' % (selection))
        _self.show_as('nonbonded', '(%s) and inorganic' % (selection))
    else:
        _self.show_as('cartoon', selection)
        _self.show_as('sticks', '(%s) and organic' % (selection))
        _self.show_as('spheres', '(%s) and inorganic' % (selection))
        _self.show_as('nonbonded', '(%s) and solvent' % (selection))


@cmd.extend
def show_ptm(selection="*", *, rep="sticks", cys=True, _self=cmd):
    """
DESCRIPTION

    Show post-translational modifications.
    """
    seleexpr = (f"({selection}) and "
                "(not (hydro and bound_to elem C)) and "
                "(byres ("
                "(bound_to not polymer) or "
                "(polymer and bound_to hetatm) or "
                "(name SG and bound_to name SG)"
                "))")
    _self.show(rep, seleexpr)

    if int(cys):
        seleexpr = f"({selection}) and name SG and bound_to name SG"
        _self.show("spheres", seleexpr)
        _self.color("yellow", seleexpr)


def get_color_family(color):
    '''
DESCRIPTION

    API only. Get the colors list from a color submenu (like all "Reds").
    '''
    from pymol import menu
    colors = []
    try:
        for (_, col) in dict(menu.all_colors_list)[color]:
            colors.append(col)
    except KeyError:
        raise CmdException(repr(color)) from None
    return colors


def cbm(selection='all', first_color=2, *, _self=cmd):
    '''
DESCRIPTION

    Color by molecule

USAGE

    cbm [ selection [, first_color ]]
    '''
    import itertools
    if isinstance(first_color, str) and first_color.endswith('s'):
        colors = get_color_family(first_color)
        col_it = itertools.cycle(colors)
    else:
        col_it = itertools.count(int(first_color))
    for model in _self.get_object_list('(' + selection + ')'):
        col = next(col_it)
        if selection in ['*', 'all']:
            # explicit model coloring for cmd.get_object_color_index to work
            _self.color(col, model)
        else:
            _self.color(col, '%s and (%s)' % (model, selection))


def cbs(selection='all', first_color=2, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Color by segment
    '''
    import itertools
    if isinstance(first_color, str) and first_color.endswith('s'):
        colors = get_color_family(first_color)
        col_it = itertools.cycle(colors)
    else:
        col_it = itertools.count(int(first_color))
    segi_colors = dict()

    def callback(segi):
        if segi not in segi_colors:
            segi_colors[segi] = next(col_it)
        return segi_colors[segi]
    _self.alter(selection, 'color = callback(segi)', space=locals())
    _self.rebuild()


expression_sc = cmd.Shortcut([
    'count',
    'resi',
    'b',
    'q',
    'pc',
])


def spectrumany(expression, color_list, selection='(all)', minimum=None, maximum=None, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Define a color spectrum with as many color-stops as you like (at least 2).

    Note: Since PyMOL 1.6, the regular `spectrum` command also supports
    arbitrary color lists, so this command is obsolete.

USAGE

    spectrumany expression, color_list [, selection [, minimum [, maximum ]]]

ARGUMENTS

    expression = count, resi, b, q, or pc: respectively, atom count, residue
    index, temperature factor, occupancy, or partial charge {default: count}

    color_list = string: Space separated list of colors

    ... all other arguments like with `spectrum` command

EXAMPLE

    spectrumany count, forest green yellow white
    spectrumany b, red yellow white, (polymer), maximum=100.0

SEE ALSO

    spectrum
    '''
    quiet = int(quiet)
    colors = color_list.split()
    if len(colors) < 2:
        print('failed! please provide at least 2 colors')
        return

    colvec = [_self.get_color_tuple(i) for i in colors]
    parts = len(colvec) - 1

    expression = {
        'pc': 'partial_charge',
        'fc': 'formal_charge',
        'count': 'index',
    }.get(expression, expression)
    minmax_expr = {'resi': 'resv'}.get(expression, expression)
    discrete_expr = ['index', 'resi']

    if _self.count_atoms(selection) == 0:
        print('empty selection')
        return

    if None in [minimum, maximum]:
        e_list = list()
        _self.iterate(selection, 'e_list.append(%s)' % (minmax_expr), space=locals())
        if minimum is None:
            minimum = min(e_list)
        if maximum is None:
            maximum = max(e_list)
    minimum, maximum = float(minimum), float(maximum)
    if not quiet:
        print(' Spectrum: range (%.5f to %.5f)' % (minimum, maximum))

    if maximum == minimum:
        print('no spectrum possible, only equal %s values' % (expression))
        return

    if expression in discrete_expr:
        val_range = int(maximum - minimum + 1)
    else:
        val_range = maximum - minimum
        _self.color(colors[0], selection)

    steps = 60 // parts
    steps_total = steps * parts

    val_start = minimum
    for p in range(parts):
        for i in range(steps):
            ii = float(i) / steps
            col_list = [colvec[p + 1][j] * ii + colvec[p][j] * (1.0 - ii) for j in range(3)]
            col_name = '0x%02x%02x%02x' % tuple(int(0xFF * v) for v in col_list)
            val_end = val_range * (i + 1 + p * steps) // steps_total + minimum
            if expression in discrete_expr:
                _self.color(col_name, '(%s) and %s %d-%d' % (selection, expression, val_start, val_end))
            else:
                _self.color(col_name, '(%s) and %s > %f' % (selection, expression, val_start))
            val_start = val_end


def spectrum_states(selection='all', representations='cartoon ribbon',
        color_list='blue cyan green yellow orange red',
        first=1, last=0, quiet=1, *, _self=cmd):
    '''
DESCRIPTION

    Color each state in a multi-state object different.

    (c) 2011 Takanori Nakane and Thomas Holder

USAGE

    spectrum_states [ selection [, representations [, color_list [, first [, last ]]]]]

ARGUMENTS

    selection = string: object names (works with complete objects only)
    {default: all}

    representations = string: space separated list of representations
    {default: cartoon ribbon}

    color_list = string: space separated list of colors {default: blue cyan
    green yellow orange red}

SEE ALSO

    spectrum, spectrumany
    '''
    from math import floor, ceil

    first, last, quiet = int(first), int(last), int(quiet)
    colors = color_list.split()
    if len(colors) < 2:
        raise CmdException('please provide at least 2 colors')

    colvec = [_self.get_color_tuple(i) for i in colors]

    # filter for valid <repr>_color settings
    settings = []
    for r in representations.split():
        if r[-1] == 's':
            r = r[:-1]
        s = r + '_color'
        if s in _self.setting.name_list:
            settings.append(s)
        elif not quiet:
            print(' Warning: no such setting: ' + repr(s))

    # object names only
    selection = ' '.join(_self.get_object_list('(' + selection + ')'))
    if _self.count_atoms(selection) == 0:
        raise CmdException('empty selection')

    if last < 1:
        last = _self.count_states(selection)

    val_range = int(last - first + 1)
    if val_range < 2:
        raise CmdException('no spectrum possible, need more than 1 state')

    for i in range(val_range):
        p = float(i) / (val_range - 1) * (len(colvec) - 1)
        p0, p1 = int(floor(p)), int(ceil(p))
        ii = (p - p0)
        col_list = [colvec[p1][j] * ii + colvec[p0][j] * (1.0 - ii) for j in range(3)]
        col_name = '0x%02x%02x%02x' % tuple(int(0xFF * v) for v in col_list)
        for s in settings:
            _self.set(s, col_name, selection, state=i + first)


class scene_preserve(object):
    '''
DESCRIPTION

    API only. Context manager to restore the current scene on exit.
    '''

    def __init__(self, *, _self=cmd, **kwargs):
        self._self = _self
        self.kwargs = kwargs

    def __enter__(self):
        import random
        self.name = 'tmp_%d' % (random.randint(0, 100_000_000))
        self._self.scene(self.name, 'store', **self.kwargs)

    def __exit__(self, type, value, traceback):
        self._self.scene(self.name, 'recall')
        self._self.scene(self.name, 'delete')


def goodsell_lighting(*, _self=cmd):
    """
DESCRIPTION

    Mostly flat lighting, similar to the famous illustrations by David
    Goodsell. Works well with spheres representation.
    """
    _self.set("ambient", 0.6)
    _self.set("direct", 0.6)
    _self.set("reflect", 0.1)
    _self.set("shininess", 60.0)
    _self.set("ray_shadow", 'off')
    _self.set("specular_intensity", 0.0)


@cmd.extend
def grid_labels(expression: str = "model",
                selection: str = "enabled",
                *,
                subsele: str = "all",
                _self=cmd):
    """
DESCRIPTION

    Create screen stabilized labels per grid cell
    """
    if expression == "title":
        expression = "cmd.get_title(model, -1)"

    obj2group = {
        obj: group
        for group in cmd.get_names_of_type("object:group")
        for obj in cmd.get_object_list(group)
    }

    slots_done: set = set()  # set[int]
    names_done: set = set()  # set[str]

    def done(set_: set, name) -> bool:
        if name in set_:
            return True
        set_.add(name)
        return False

    tmpsele = _self.get_unused_name("_tmpsele")

    for i, model in enumerate(_self.get_object_list(selection), 1):
        name = obj2group.get(model, model)
        if done(names_done, name):
            continue

        slot = _self.get_setting_int("grid_slot", name)
        if done(slots_done, i if slot == -1 else slot):
            continue

        _self.select(tmpsele, f"first (%{model} and ({subsele}))")
        try:
            _self.label(tmpsele, expression)
            _self.set("label_screen_point", (-1, 1, 0), tmpsele)
            _self.set("label_position", (1, -1, 0), tmpsele)
            _self.set("label_relative_mode", 1, tmpsele)
        finally:
            _self.delete(tmpsele)


# commands
cmd.alias('z', 'zoom visible')
cmd.alias('x', 'nice')
cmd.extend('nice', nice)
cmd.extend('cbm', cbm)
cmd.extend('cbs', cbs)
cmd.extend('spectrumany', spectrumany)
cmd.extend('spectrum_states', spectrum_states)
cmd.extend('goodsell_lighting', goodsell_lighting)

# tab-completion of arguments
cmd.auto_arg[0]['show_ptm'] = cmd.auto_arg[2]['spectrum']
cmd.auto_arg[0]['spectrumany'] = [expression_sc, 'expression', ', ']
cmd.auto_arg[1]['spectrumany'] = [cmd.auto_arg[0]['color'][0], 'color', ' ']
cmd.auto_arg[2]['spectrumany'] = cmd.auto_arg[2]['spectrum']
cmd.auto_arg[0]['spectrum_states'] = cmd.auto_arg[0]['disable']
cmd.auto_arg[1]['spectrum_states'] = [cmd.auto_arg[0]['show'][0], 'representation', ' ']
cmd.auto_arg[2]['spectrum_states'] = cmd.auto_arg[1]['spectrumany']

# vi:expandtab:smarttab
