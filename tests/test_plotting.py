import psico.plotting
import matplotlib.pyplot
from matplotlib.testing.decorators import image_comparison
from pathlib import Path
from pymol import cmd

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_get_model_color():
    cmd.reinitialize()
    cmd.pseudoatom("m1")
    cmd.color("0x123456", "m1")
    assert psico.plotting.get_model_color("m1") == "#123456"


# def rms_plot(selection='guide', ref1=None, ref2=None, state1=1, state2=-1,
# def area_plot(selection='all', filename=None, *, quiet=1, _self=cmd):
# def pca_plot(aln_object, ref='all', state=0, maxlabels=20, size=20, invert='',
# def iterate_plot(selection, expr_y, expr_x=None, scatter=0, filename=None,


@image_comparison(baseline_images=['contact_map_plot-2x19-frag-mse.png'],
                  remove_text=True,
                  style='mpl20')
def test_contact_map_plot():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "2x19-frag-mse.pdb", "m1")
    matplotlib.pyplot.rcParams['figure.figsize'] = [0.64, 0.48]
    psico.plotting.contact_map_plot()
