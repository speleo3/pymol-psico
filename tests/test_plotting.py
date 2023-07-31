import psico.plotting
import functools
import matplotlib.pyplot
from matplotlib.testing.decorators import image_comparison
from pathlib import Path
from pymol import cmd

DATA_PATH = Path(__file__).resolve().parent / 'data'


def my_image_comparison(key: str = "default", figsize=(0.64, 0.48)):

    def decorator(func):
        assert func.__name__.startswith("test_")

        @image_comparison([f'{func.__name__[5:]}-{key}.png'],
                          remove_text=True,
                          style='mpl20')
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            cmd.reinitialize()
            matplotlib.pyplot.rcParams['figure.figsize'] = figsize
            return func(*args, **kwargs)

        return wrapper

    return decorator


def test_get_model_color():
    cmd.reinitialize()
    cmd.pseudoatom("m1")
    cmd.color("0x123456", "m1")
    assert psico.plotting.get_model_color("m1") == "#123456"


@my_image_comparison()
def test_rms_plot():
    cmd.load(DATA_PATH / "1nmr-frag-nohydro.pdb", "m1")
    cmd.color("blue")
    psico.plotting.rms_plot()


@my_image_comparison()
def test_area_plot():
    cmd.load(DATA_PATH / "1nmr-frag-nohydro.pdb", "m1")
    cmd.color("red")
    psico.plotting.area_plot()


def test_pca_plot(tmp_path):
    # PCA plot can be randomly inverted (same PCA), so don't use image comparison.
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1nmr-frag-nohydro.pdb", "m1")
    cmd.split_states("m1")
    cmd.extra_fit("m1_*", object="aln")
    path = tmp_path / "plot.png"
    psico.plotting.pca_plot("aln", filename=path)
    assert path.exists()


@my_image_comparison()
def test_iterate_plot():
    cmd.load(DATA_PATH / "2x19-frag-mse.pdb", "m1")
    cmd.color("magenta")
    psico.plotting.iterate_plot("guide", "b", "resv")


@my_image_comparison('2x19-frag-mse')
def test_contact_map_plot():
    cmd.load(DATA_PATH / "2x19-frag-mse.pdb", "m1")
    psico.plotting.contact_map_plot()
