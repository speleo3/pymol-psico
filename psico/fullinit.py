'''
Convenience module for single line psico initialization

(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

if __name__.endswith('.fullinit'):
    from . import init
    from . import guitweak  # noqa: F401
    from . import pymol_version

    init(pymol_version < 2.5, False, 1)
else:
    import os, imp

    try:
        _script_path = __script__  # noqa: F821 Undefined name
    except NameError:
        raise ImportError('invalid invocation of psico.fullinit')

    imp.load_module('psico', None, os.path.dirname(_script_path),
            ('', '', imp.PKG_DIRECTORY))

    import psico.fullinit  # noqa: F401

# vi: ts=4:sw=4:smarttab:expandtab
