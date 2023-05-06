'''
Collection of PyMOL scripts

A reference document with all commands can be generated with:
PyMOL> psico.helping.write_html_ref('psico-commands.html')

(c) 2010-2012 Thomas Holder <speleo3@users.sourceforge.net>
              Steffen Schmidt <steffen.schmidt@tuebingen.mpg.de>
              Max Planck Institute for Developmental Biology

License: BSD-2-Clause
'''

__version__ = '4.2'

from psico.versioning import make_version_int_tuple

try:
    from pymol import cmd
    pymol_version_str, pymol_version = cmd.get_version()[:2]
except Exception as ex:
    print(ex)
    pymol_version_tuple = (2, 0)
    pymol_version = pymol_version_tuple[0] + pymol_version_tuple[1] * 0.1
    pymol_version_str = "%d.%d" % pymol_version_tuple
else:
    pymol_version_tuple = make_version_int_tuple(pymol_version_str)

__all__ = [
    'aaindex',
    'aggrescanning',
    'conservation',
    'creating',
    'editing',
    'electrostatics',
    'exporting',
    'fasta',
    'fitting',
    'geometry',
    'orientation',
    'helping',
    'importing',
    'minimizing',
    'modelling',
    'moving',
    'msms',
    'nma',
    'plotting',
    'querying',
    'selecting',
    'setting',
    'snp',
    'viewing',
    'wizards',
    'xtal',
    'load_mtz_cctbx',
]


def make_global():
    '''
    "psico" might be installed as submodule of something else (PyMOL plugin).
    Invoke this function if you want to do "from psico import ...".
    '''
    import sys
    if sys.modules.get('psico') != sys.modules[__name__]:
        sys.modules['psico'] = sys.modules[__name__]


def init(save=0, fetch=0, pymolapi=0):
    '''
DESCRIPTION

    Imports all psico submodules and puts "psico" into the pymol namespace
    for GUI menus. Also enables "help psico" in the PyMOL command line.

ARGUMENTS

    save = bool: Overload "save" command with psico.exporting.save (writes
    secondary structure and crystal records to PDB header)

    fetch = bool: Overload "fetch" command with psico.importing.fetch (can
    fetch from local mirror, knows SCOP and CATH identifiers)

    pymolapi = 0/1/2: Add all psico functions to PyMOL API (pymol.cmd)
    '''
    import pymol
    from pymol import cmd

    # init all submodules
    psico = __import__(__name__, fromlist=__all__)

    # pymol namespace
    if not hasattr(pymol, 'psico'):
        pymol.psico = psico

    # pymol help
    if 'psico' not in cmd.help_only:
        cmd.help_only['psico'] = [psico]
        cmd.help_sc.append('psico')

    if save:
        cmd.extend('save', psico.exporting.save)

    if fetch:
        cmd.extend('fetch', psico.importing.fetch)

    if pymolapi:
        init_cmd(pymolapi == 2)

# PyMOL Plugin hook


def __init_plugin__(self=None):
    init(1, 0, 0)
    make_global()


def init_cmd(force=0):
    '''
    Adds all psico functions to PyMOL API (pymol.cmd)

    If force is True, overwrite existing names.
    '''
    from pymol import cmd

    for name, value in cmd.keyword.items():
        function = value[0]
        if function.__module__.startswith(__name__):
            if force or not hasattr(cmd, function.__name__):
                setattr(cmd, function.__name__, function)


# See also http://pymolwiki.org/index.php/Aa_codes
one_letter = {
    'PAQ': 'Y', 'AGM': 'R', 'ILE': 'I', 'PR3': 'C', 'GLN': 'Q', 'DVA': 'V',
    'CCS': 'C', 'ACL': 'R', 'GLX': 'Z', 'GLY': 'G', 'GLZ': 'G', 'DTH': 'T',
    'OAS': 'S', 'C6C': 'C', 'NEM': 'H', 'DLY': 'K', 'MIS': 'S', 'SMC': 'C',
    'GLU': 'E', 'NEP': 'H', 'BCS': 'C', 'ASQ': 'D', 'ASP': 'D', 'SCY': 'C',
    'SER': 'S', 'LYS': 'K', 'SAC': 'S', 'PRO': 'P', 'ASX': 'B', 'DGN': 'Q',
    'DGL': 'E', 'MHS': 'H', 'ASB': 'D', 'ASA': 'D', 'NLE': 'L', 'DCY': 'C',
    'ASK': 'D', 'GGL': 'E', 'STY': 'Y', 'SEL': 'S', 'CGU': 'E', 'ASN': 'N',
    'ASL': 'D', 'LTR': 'W', 'DAR': 'R', 'VAL': 'V', 'CHG': 'A', 'TPO': 'T',
    'CLE': 'L', 'GMA': 'E', 'HAC': 'A', 'AYA': 'A', 'THR': 'T', 'TIH': 'A',
    'SVA': 'S', 'MVA': 'V', 'SAR': 'G', 'LYZ': 'K', 'BNN': 'A', '5HP': 'E',
    'IIL': 'I', 'SHR': 'K', 'HAR': 'R', 'FME': 'M', 'PYX': 'C', 'ALO': 'T',
    'PHI': 'F', 'ALM': 'A', 'PHL': 'F', 'MEN': 'N', 'TPQ': 'A', 'GSC': 'G',
    'PHE': 'F', 'ALA': 'A', 'MAA': 'A', 'MET': 'M', 'UNK': 'X', 'LEU': 'L',
    'ALY': 'K', 'SET': 'S', 'GL3': 'G', 'TRG': 'K', 'CXM': 'M', 'TYR': 'Y',
    'SCS': 'C', 'DIL': 'I', 'TYQ': 'Y', '3AH': 'H', 'DPR': 'P', 'PRR': 'A',
    'CME': 'C', 'IYR': 'Y', 'CY1': 'C', 'TYY': 'Y', 'HYP': 'P', 'DTY': 'Y',
    '2AS': 'D', 'DTR': 'W', 'FLA': 'A', 'DPN': 'F', 'DIV': 'V', 'PCA': 'E',
    'MSE': 'M', 'MSA': 'G', 'AIB': 'A', 'CYS': 'C', 'NLP': 'L', 'CYQ': 'C',
    'HIS': 'H', 'DLE': 'L', 'CEA': 'C', 'DAL': 'A', 'LLP': 'K', 'DAH': 'F',
    'HMR': 'R', 'TRO': 'W', 'HIC': 'H', 'CYG': 'C', 'BMT': 'T', 'DAS': 'D',
    'TYB': 'Y', 'BUC': 'C', 'PEC': 'C', 'BUG': 'L', 'CYM': 'C', 'NLN': 'L',
    'CY3': 'C', 'HIP': 'H', 'CSO': 'C', 'TPL': 'W', 'LYM': 'K', 'DHI': 'H',
    'MLE': 'L', 'CSD': 'A', 'HPQ': 'F', 'MPQ': 'G', 'LLY': 'K', 'DHA': 'A',
    'DSN': 'S', 'SOC': 'C', 'CSX': 'C', 'OMT': 'M', 'DSP': 'D', 'PTR': 'Y',
    'TRP': 'W', 'CSW': 'C', 'EFC': 'C', 'CSP': 'C', 'CSS': 'C', 'SCH': 'C',
    'OCS': 'C', 'NMC': 'G', 'SEP': 'S', 'BHD': 'D', 'KCX': 'K', 'SHC': 'C',
    'C5C': 'C', 'HTR': 'W', 'ARG': 'R', 'TYS': 'Y', 'ARM': 'R', 'DNP': 'A',
}
three_letter = {
    'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE',
    'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN',
    'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP',
    'V': 'VAL', 'Y': 'TYR',
}

# matplotlib_fix preferences
matplotlib_fix_prefs = {
    'verbose': True,
    'force_tkagg': False,
    'force_show': False,
    'tkagg_overload': True,
}


def which(*names, **kw):
    '''
    Return full path to executable or empty string if not found in PATH.
    '''
    import os
    if 'PATHEXT' in os.environ:
        pathext = [''] + os.environ['PATHEXT'].split(';')
        names = [n + ext for n in names for ext in pathext]
    path = kw.get('path') or os.environ['PATH'].split(os.pathsep)
    for n in names:
        for p in path:
            full = os.path.join(p, n)
            if os.path.isfile(full):
                return full
    return ''

# vi: ts=4:sw=4:smarttab:expandtab
