import os
import sys
import textwrap
from unittest.mock import NonCallableMock, NonCallableMagicMock

# -- Path setup --------------------------------------------------------------
# Add project root to sys.path so Sphinx can import the package
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)

# -- Minimal mocks for heavy/optional dependencies ---------------------------
# Mock the PyMOL environment and other heavy optional imports so autodoc works
# on Read the Docs without requiring binary dependencies.
mock_module_names = [
    'pymol', 'pymol.movie', 'pymol.wizard',
    'chempy',
    'numpy',
    'scipy',
    'Bio', 'Bio.Seq', 'Bio.Align', 'Bio.PDB',
    'csb',
    'prody',
    'rdkit', 'rdkit.Chem', 'rdkit.Chem.AllChem', 'rdkit.DataStructs',
    'openbabel',
    'mdtraj',
    'matplotlib', 'matplotlib.pyplot', 'matplotlib.backends', 'matplotlib.backends.backend_agg',
    'epam', 'epam.indigo',
]

for name in mock_module_names:
    if name not in sys.modules:
        sys.modules[name] = NonCallableMock(name=name)


class CmdException(Exception):
    pass


class _CmdMock(NonCallableMagicMock):

    def extend(self, name, function=None):
        return function or name

    def extendaa(self, *args, **kwargs):
        return self.extend

    def get_version(self):
        return ("2.5.0", 2.5, 0)

    def __repr__(self) -> str:
        # function signatures render _self=...
        return '...'


pymol_mock = sys.modules["pymol"]
pymol_mock.cmd = _CmdMock()  # type: ignore[attr-defined]
pymol_mock.CmdException = CmdException  # type: ignore[attr-defined]


def process_pymol_docstring(app, what, name: str, obj, options,
                            lines: list[str]):
    """
    Process PyMOL-style docstring to strip the 'DESCRIPTION' header, dedent
    the description body, and apply basic reStructuredText markup.
    """
    for i, line in enumerate(lines):
        if line not in ('', 'DESCRIPTION'):
            lines[:i] = []
            break

    j = 0
    for i, line in enumerate(lines, 1):
        if line:
            if not line[0].isspace():
                break
            j = i

    lines[:j] = textwrap.dedent("\n".join(lines[:j])).splitlines()

    for i in range(len(lines), j, -1):
        line = lines[i - 1]
        if line.isupper() and not line.startswith(" "):
            lines[i - 1] = line.title()
            lines.insert(i, "~" * len(line))


# -- Project information -----------------------------------------------------
project = 'psico'
copyright = 'BSD-2-Clause'
author = 'Thomas Holder and contributors'

# The full version, including alpha/beta/rc tags
try:
    from psico import __version__ as release
except Exception:
    release = 'unknown'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

autosummary_generate = True
napoleon_google_docstring = True
napoleon_numpy_docstring = True

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': False,
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']

# -- Intersphinx -------------------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

# -- Setup function for autodoc processor -----------------------------------
def setup(app):
    app.connect('autodoc-process-docstring', process_pymol_docstring)
