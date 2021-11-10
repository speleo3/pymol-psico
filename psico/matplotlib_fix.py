'''
PyMOL matplotlib integration

matplotlib uses Tk, which conflicts with the PyMOL Tk mainloop. This module
overloads the TkAgg backend to use PyMOLs Tk instance as master widget.

This is likely to be fragile against matplotlib upgrades.

(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

def overload():
    import matplotlib
    from . import matplotlib_fix_prefs

    if matplotlib_fix_prefs['verbose']:
        print(' matplotlib_fix: If you have issues with matplotlib and PyMOL, check')
        print(' matplotlib_fix: the values of psico.matplotlib_fix_prefs')

    if matplotlib_fix_prefs['force_tkagg']:
        matplotlib.use('TkAgg')

    if matplotlib.get_backend() != 'TkAgg':
        return

    if not matplotlib_fix_prefs['tkagg_overload']:
        return

    from matplotlib.backends import backend_tkagg

    def _new_figure_manager(num, *args, **kwargs):
        import pymol

        if pymol._ext_gui is None:
            return new_figure_manager(num, *args, **kwargs)

        backend_tkagg.show._needmain = False

        import tkinter as Tk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureManagerTkAgg

        FigureClass = kwargs.pop('FigureClass', Figure)
        figure = FigureClass(*args, **kwargs)

        window = Tk.Toplevel(master=pymol._ext_gui.root)

        canvas = FigureCanvasTkAgg(figure, master=window)
        figManager = FigureManagerTkAgg(canvas, num, window)
        if matplotlib.is_interactive():
            figManager.show()
        return figManager

    new_figure_manager = backend_tkagg.new_figure_manager
    backend_tkagg.new_figure_manager = _new_figure_manager

try:
    overload()
except:
    print('matplotlib_fix failed')

# vi:expandtab:smarttab
