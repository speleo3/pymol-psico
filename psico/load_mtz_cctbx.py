'''
(c) Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
'''

import os
import sys


def mtz2ccp4maps(filename, prefix, amplitudes, phases):
    '''
    Creates a temporary directory and dumps all maps from the given MTZ file
    into this directory as CCP4 maps files. Returns the path of the temporary
    directory.
    '''
    import tempfile
    from iotbx.reflection_file_reader import any_reflection_file

    hkl_in = any_reflection_file(file_name=filename)

    temp_dir = tempfile.mkdtemp()

    for array in hkl_in.as_miller_arrays():
        if not array.is_complex_array():
            continue

        labels = array.info().labels
        if amplitudes and [amplitudes, phases] != labels:
            continue

        fft_map = array.fft_map(resolution_factor=0.25).apply_sigma_scaling()
        map_filename = os.path.join(temp_dir, prefix + '_' + '_'.join(labels))
        if hasattr(fft_map, 'as_ccp4_map'):
            fft_map.as_ccp4_map(file_name=map_filename + '.ccp4')
        else:
            fft_map.as_xplor_map(file_name=map_filename + '.xplor')

    return temp_dir


if __name__ == '__main__':
    # Standalone script:
    # print the name of the temporary directory to standard output

    print(mtz2ccp4maps(*sys.argv[1:]))
    sys.exit(0)

try:
    # running this script in the "pymol" namespace
    this_file = __script__  # type: ignore[name-defined]
except NameError:
    this_file = __file__

import pymol  # noqa: E402


def load_mtz_cctbx(filename, prefix='', amplitudes='', phases='', quiet=1,
        _self=pymol.cmd):
    '''
DESCRIPTION

    Load maps from an MTZ file, using iotbx (via "cctbx.python").

    Map objects will be named: <prefix>_<amplitudes>_<phases>

ARGUMENTS

    filename = str: path to mtz file

    prefix = str: object name prefix for new map objects

    amplitudes = str: amplitudes column label (optional). If not given,
    load all maps. If given, the 'phases' argument is required as well.

    phases = str: phases column label (required if and only if
    amplitudes is given as well)
    '''
    import subprocess
    import shutil

    if not prefix:
        prefix = os.path.basename(filename).rpartition('.')[0]

    args = [filename, prefix, amplitudes, phases]

    try:
        # try with "cctbx.python"
        process = subprocess.Popen(['cctbx.python', this_file] + args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            raise pymol.CmdException(stderr)
        outdir = stdout.strip()
        if not isinstance(outdir, str):
            outdir = outdir.decode()
    except OSError:
        try:
            # try inside this Python interpreter
            outdir = mtz2ccp4maps(*args)
        except ImportError as ex:
            raise pymol.CmdException("can't import iotbx and can't run cctbx.python") from ex

    # normalization is done by apply_sigma_scaling()
    normalize = _self.get_setting_int('normalize_ccp4_maps')
    if normalize:
        _self.set('normalize_ccp4_maps', 0)

    for mapfilename in os.listdir(outdir):
        _self.load(os.path.join(outdir, mapfilename), quiet=quiet)

    if normalize:
        _self.set('normalize_ccp4_maps', normalize)

    shutil.rmtree(outdir)


pymol.cmd.extend('load_mtz_cctbx', load_mtz_cctbx)
