import psico.fitting
import psico.querying
import os
from pymol import cmd
from pytest import approx

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
FILENAME_MULTISTATE = os.path.join(DATA_DIR, '1nmr-frag-nohydro.pdb')


def test_intra_xfit():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.fitting.intra_xfit("all", load_b=1, cycles=10)
    b_list = psico.querying.iterate_to_list("guide", "b")
    b_list_ref = [
        -3.932, -3.502, -3.806, -3.425, -4.862, -5.246, -4.358, -3.846, -4.348,
        -5.397, -4.353, -4.296, -4.037, -4.310, -4.874, -5.439, -4.552, -3.804,
        -3.987, -3.468, -3.019, -3.061, -1.961, -1.850, -1.721, -1.445, -1.088,
        -0.109, 2.399, 1.314, 2.389
    ]
    assert b_list == approx(b_list_ref, abs=1e-3)


def test_intra_center():
    cmd.reinitialize()
    cmd.load(FILENAME_MULTISTATE)
    psico.fitting.intra_center("resi 80")
    assert cmd.intra_rms_cur("*") == approx(
        [-1.0, 1.425962, 8.287644, 7.676795])


# def alignwithanymethod(mobile, target, methods=None, async_=1, quiet=1, **kwargs):
# def tmalign(mobile, target, mobile_state=1, target_state=1, args='',
# def dyndom_parse_info(filename, selection='(all)', quiet=0):
# def dyndom(mobile, target, window=5, domain=20, ratio=1.0, exe='', transform=1,
# def gdt_ts(mobile, target, cutoffs='1 2 4 8', quiet=1):
# def get_rmsd_func():
# def matchmaker(mobile, target, match):
# def local_rms(mobile, target, window=20, mobile_state=1, target_state=1,
# def extra_fit(selection='(all)', reference=None, method='align', zoom=1,
# def theseus(mobile, target, match='align', cov=0, cycles=200,
# def intra_theseus(selection, state=1, cov=0, cycles=200,
# def prosmart(mobile, target, mobile_state=1, target_state=1,
# def xfit(mobile, target, mobile_state=-1, target_state=-1, load_b=0,
# def promix(mobile, target, K=0, prefix=None, mobile_state=-1, target_state=-1,
# def intra_promix(selection, K=0, prefix=None, conformers=0, guide=1,
# def intra_boxfit(selection="polymer", center=[0.5, 0.5, 0.5], _self=cmd):
