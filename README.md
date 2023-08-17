[![CI](https://github.com/speleo3/pymol-psico/workflows/CI/badge.svg)](https://github.com/speleo3/pymol-psico/actions)
[![codecov](https://codecov.io/gh/speleo3/pymol-psico/branch/master/graph/badge.svg)](https://codecov.io/gh/speleo3/pymol-psico)
[![Anaconda.org](https://anaconda.org/speleo3/pymol-psico/badges/version.svg)](https://anaconda.org/speleo3/pymol-psico)

# Pymol ScrIpt COllection (PSICO)

*psico* is a python module which extends PyMOL, a molecular visualization
tool (https://pymol.org).

## Documentation

https://pymolwiki.org/index.php/psico

A reference document with all commands can be generated with:

```python
psico.helping.write_html_ref('psico-commands.html')
```

## Installation

Option 1: Install from source. In your activated (conda or venv) environment, run:

```sh
pip install .
```

Option 2: Use a [conda-compatible](https://mamba.readthedocs.io/) package manager. In your activated conda environment, run:

```sh
micromamba install speleo3::pymol-psico
```

To activate all *psico* commands in PyMOL, add this to your `~/.pymolrc.py`
file (or just enter into the PyMOL command line):

```python
import psico.fullinit
```

## Dependencies

* PyMOL 2.3+ (https://pymol.org/)
* Python 3.7+

Some functions in *psico* require the following python modules:

* numpy
* Bio (biopython)
* csb (https://github.com/csb-toolbox/CSB)
* modeller
* prody
* rdkit (http://www.rdkit.org/)
* indigo (https://github.com/epam/Indigo)

*psico* has several wrappers for external tools. For full support, you need
these binaries:

* dssp (https://swift.cmbi.umcn.nl/gv/dssp/)
* stride (http://webclu.bio.wzw.tum.de/stride/)
* TMalign (https://zhanggroup.org/TM-align/)
* TMscore
* MMalign
* theseus (http://www.theseus3d.org)
* needle (http://emboss.sourceforge.net)
* DynDom (http://dyndom.cmp.uea.ac.uk/dyndom/)
* qdelaunay (http://www.qhull.org)
* mencoder (http://www.mplayerhq.hu)
* pdbmat (http://ecole.modelisation.free.fr/modes.html)
* diagrtb
* apbs (http://www.poissonboltzmann.org/)
* pdb2pqr
* prosmart (http://www2.mrc-lmb.cam.ac.uk/groups/murshudov/)
* p_sstruc3
* msms (https://ccsb.scripps.edu/msms/)
