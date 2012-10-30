'''
Plotting with matplotlib

(c) 2011-2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd, CmdException

def _showfigure(fig, filename, quiet):
    '''
DESCRIPTION

    Helper function for plot commands.
    '''
    if filename is None:
        try:
            fig.show()
        except AttributeError:
            print ' Warning: matplotlib.pyplot.figure.show failed'
            filename, quiet = 'psico-plot-failsafe.pdf', False
        else:
            from . import matplotlib_fix_prefs
            if matplotlib_fix_prefs['force_show']:
                import matplotlib
                matplotlib.pyplot.show()
            return
    try:
        fig.savefig(filename)
        if not quiet:
            print ' Plot written to', filename
    except ValueError, e:
        print ' Error:', e

def get_model_color(model):
    '''
DESCRIPTION

    API only. Get model color as #xxxxxx to be used with matplotlib.
    '''
    from .querying import get_color
    return get_color(model, 0, 2)

def rms_plot(selection='guide', ref1=None, ref2=None, state1=1, state2=-1,
        match='align', cur=0, maxlabels=20, size=20, alpha=0.75,
        filename=None, quiet=1):
    '''
DESCRIPTION

    Scatterplot with RMSD for all models/states in selection against ref1 and
    ref2.

ARGUMENTS

    selection = string: atom selection {default: guide}

    ref1 = string: atom selection {default: first model in selection}
    ref2 = string: atom selection {default: last model in selection}

    state1 = int: object state for ref1
    state2 = int: object state for ref2

    match = string: in, like, align, none or the name of an alignment object
    (see "local_rms" help for details) {default: align}

    cur = 0/1: if 1, use rms_cur instead of rms (no fitting) {default: 0}
    '''
    from .fitting import matchmaker
    from . import matplotlib_fix
    from matplotlib.pyplot import figure

    state1, state2, cur, quiet = int(state1), int(state2), int(cur), int(quiet)
    maxlabels, size, alpha = int(maxlabels), float(size), float(alpha)

    models = cmd.get_object_list('(' + selection + ')')
    if ref1 is None:
        ref1 = models[0]
    if ref2 is None:
        ref2 = models[-1]
    if state2 == -1:
        if ref1 == ref2:
            state2 = cmd.count_states(ref2)
        else:
            state2 = 1

    sele_pairs = {}
    tmp_names_all = []
    for model in models:
        for ref in [ref1, ref2]:
            key = (model, ref)
            if key in sele_pairs:
                continue
            if model == ref:
                mobile = target = '(%s) and (%s)' % (model, selection)
            else:
                mobile, target, tmp_names = matchmaker('(%s) and (%s)' % \
                        (model, selection), ref, match)
                tmp_names_all.extend(tmp_names)
            sele_pairs[key] = (mobile, target)

    x_list = []
    y_list = []
    colors = []
    text_list = []

    _rms = cmd.rms_cur if cur else cmd.rms
    def rms(model, ref, state, ref_state):
        mobile, target = sele_pairs[model, ref]
        return _rms(mobile, target, state, ref_state, matchmaker=4)

    for model in models:
        for state in range(1, cmd.count_states(model)+1):
            rms1 = rms(model, ref1, state, state1)
            rms2 = rms(model, ref2, state, state2)
            x_list.append(rms1)
            y_list.append(rms2)
            colors.append(get_model_color(model))
            text_list.append('%s(%d)' % (model, state))

    for name in tmp_names_all:
        cmd.delete(name)

    fig = figure()
    plt = fig.add_subplot(111, xlabel='RMSD to %s (state %d)' % (ref1, state1),
            ylabel='RMSD to %s (state %d)' % (ref2, state2),
            xlim=(0, 1+max(x_list)), ylim=(0, 1+max(y_list)))
    plt.scatter(x_list, y_list, float(size), colors, linewidths=0, alpha=float(alpha))

    if maxlabels < 0 or len(text_list) <= maxlabels:
        for (x, y, text) in zip(x_list, y_list, text_list):
            plt.text(x, y, text, horizontalalignment='left')

    _showfigure(fig, filename, quiet)

def pca_plot(aln_object, ref='all', state=0, maxlabels=20, size=20, invert='',
        which=(0,1), alpha=0.75, filename=None, quiet=1, load_b=0):
    '''
DESCRIPTION

    Principal Component Analysis on a set of superposed conformations, given
    by an alignment object. By default all states in all objects are
    considered. Generates a 2d-plot of the first two principal components.

USAGE

    pca_plot aln_object [, ref [, state [, maxlabels ]]]

ARGUMENTS

    aln_object = string: name of alignment object, defines the selection
    and the atom mapping between objects

    ref = string: object names for which to calculate PCA for {default: all}

    state = integer: if state=0 use all states {default: 0}

    maxlabels = integer: label dots in plot if maxlabels<0 or number of models
    not more than maxlabels {default: 20}

    size = float: size of plot points in px^2 {default: 20}

    invert = string: invert plotting axes x, y or xy {default: ''}

    which = (int,int): indices of principal components to plot {default: (0,1)}

    alpha = float: opacity of plotting points

    filename = string: if given, plot to file {default: None}

EXAMPLE

    fetch 1ake 4ake 1dvr 1ak2, async=0
    split_chains
    extra_fit (*_*) and name CA, reference=1ake_A, cycles=0, object=aln
    pca_plot aln, 1ake_* 4ake_*

    fetch 1ubq 2k39, async=0
    align 2k39, 1ubq and guide, cycles=0, object=aln2
    color blue, 1ubq
    color orange, 2k39
    pca_plot aln2, filename=pca-ubq.pdf
    '''
    from numpy import array, dot
    from numpy.linalg import svd, LinAlgError
    from . import matplotlib_fix
    from matplotlib.pyplot import figure

    state, quiet = int(state), int(quiet)
    maxlabels = int(maxlabels)
    if cmd.is_string(which):
        which = cmd.safe_list_eval(which)

    if aln_object not in cmd.get_names_of_type('object:'):
        print ' Warning: first argument should be an alignment object'

        from .fitting import extra_fit

        selection = aln_object
        aln_object = cmd.get_unused_name('aln')
        extra_fit(selection, cycles=0, transform=0, object=aln_object)

    if state == 0:
        states = range(1, cmd.count_states()+1)
    elif state < 0:
        states = [cmd.get_state()]
    else:
        states = [state]

    models = cmd.get_object_list(aln_object)
    references = set(cmd.get_object_list('(' + ref + ')')).intersection(models)
    others = set(models).difference(references)
    aln = cmd.get_raw_alignment(aln_object)

    if not quiet:
        print ' PCA References:', ', '.join(references)
        print ' PCA Others:', ', '.join(others)

    if len(references) == 0:
        print ' PCA Error: No reference objects'
        raise CmdException

    model_count = len(models)
    coords = dict((model, []) for model in models)
    aln = filter(lambda pos: len(pos) == model_count, aln)

    for state in states:
        idx2xyz = dict()
        cmd.iterate_state(state, aln_object, 'idx2xyz[model,index] = (x,y,z)',
                space={'idx2xyz': idx2xyz})

        for pos in aln:
            for idx in pos:
                if idx not in idx2xyz:
                    continue

                c = coords[idx[0]]
                if len(c) < state:
                    c.append([])
                c[-1].extend(idx2xyz[idx])

    c_iter = lambda models: ((c,model,i+1) for model in models
            for (i,c) in enumerate(coords[model]))
    X = array([i[0] for i in c_iter(references)])
    Y = array([i[0] for i in c_iter(others)])
    center = X.mean(0)
    X = X - center

    try:
        U, L, V = svd(X)
    except LinAlgError, e:
        print ' PCA Error: ', e
        raise CmdException

    if int(load_b):
        cmd.alter('byobj ' + aln_object, 'b=-0.01')
        b_dict = {}
        i = which[0]
        b_array = (V[i].reshape((-1, 3))**2).sum(1)**0.5
        for pos, b in zip(aln, b_array):
            for idx in pos:
                b_dict[idx] = b
        cmd.alter(aln_object, 'b=b_dict.get((model,index), -0.01)', space=locals())
        cmd.color('yellow', 'byobj ' + aln_object)
        cmd.spectrum('b', 'blue_red', aln_object + ' and b > -0.01')

    X_labels = [i[1:3] for i in c_iter(references)]
    Y_labels = [i[1:3] for i in c_iter(others)]

    x_list = []
    y_list = []
    colors = []
    text_list = []

    def plot_pc_2d(X, labels):
        pca_12 = dot(X, V.T)[:,which]
        for (x,y), (model,state) in zip(pca_12, labels):
            x_list.append(x)
            y_list.append(y)
            colors.append(get_model_color(model))
            if maxlabels < 0 or len(pca_12) <= maxlabels:
                text_list.append('%s(%d)' % (model, state))
            else:
                text_list.append(None)

    plot_pc_2d(X, X_labels)
    if len(Y) > 0:
        Y = Y - center
        plot_pc_2d(Y, Y_labels)

    if 'x' in invert:
        x_list = [-x for x in x_list]
    if 'y' in invert:
        y_list = [-y for y in y_list]

    fig = figure()
    plt = fig.add_subplot(111, xlabel='PC %d' % (which[0]+1), ylabel='PC %d' % (which[1]+1))
    plt.scatter(x_list, y_list, float(size), colors, linewidths=0, alpha=float(alpha))

    for (x, y, text) in zip(x_list, y_list, text_list):
        if text is not None:
            plt.text(x, y, text, horizontalalignment='left')

    _showfigure(fig, filename, quiet)

def iterate_plot(selection, expr_y, expr_x=None, scatter=0, filename=None,
        space=None, quiet=1):
    '''
DESCRIPTION

    Plot atomic properties.

ARGUMENTS

    selection = string: atom selection

    expr_y = string: python expression for y values
    
    expr_x = string: python expression for x values {default: None}

    scatter = 0/1: make line plot or scatter plot {default: 0, line plot}

EXAMPLE

    # C-alpha b-factors
    iterate_plot name CA, b, resv
    '''
    from . import matplotlib_fix
    from matplotlib.pyplot import figure

    scatter, quiet = int(scatter), int(quiet)
    if space is None:
        space = {'cmd': cmd, 'stored': cmd.pymol.stored}

    if cmd.is_string(selection):
        if selection.startswith('['):
            sele_list = selection[1:-1].split(',')
        else:
            sele_list = ['(%s) and (%s)' % (model, selection) for model in
                    cmd.get_object_list('(' + selection + ')')]
    else:
        sele_list = selection

    fig = figure()
    plt = fig.add_subplot(111)

    for selection in sele_list:
        space['_values'] = y_values = []
        cmd.iterate(selection, '_values.append(' + expr_y + ')', space=space)

        if expr_x is None:
            x_values = range(len(y_values))
        else:
            space['_values'] = x_values = []
            cmd.iterate(selection, '_values.append(' + expr_x + ')', space=space)

        color = get_model_color(selection)

        if scatter:
            plt.scatter(x_values, y_values, c=color)
        else:
            plt.plot(x_values, y_values, c=color)

    _showfigure(fig, filename, quiet)

# pymol commands
cmd.extend('rms_plot', rms_plot)
cmd.extend('pca_plot', pca_plot)
cmd.extend('iterate_plot', iterate_plot)

_auto_arg_aln_objects = [
    lambda: cmd.Shortcut(cmd.get_names_of_type('object:')),
    'alignment object', '']

# autocompletion
cmd.auto_arg[0].update([
    ('pca_plot', _auto_arg_aln_objects),
    ('rms_plot', cmd.auto_arg[0]['align']),
    ('iterate_plot', cmd.auto_arg[0]['iterate']),
])
cmd.auto_arg[1].update([
    ('pca_plot', cmd.auto_arg[0]['disable']),
    ('rms_plot', cmd.auto_arg[0]['align']),
])
cmd.auto_arg[2].update([
    ('rms_plot', cmd.auto_arg[1]['align']),
])

# vi: ts=4:sw=4:smarttab:expandtab
