'''
Convenience module for single line psico initialization

(c) 2012 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

if __name__.endswith('.fullinit'):
    from . import init
    from . import guitweak
    from . import pymol_version

    init(pymol_version < 2.5, False, 1)
else:
    import os, sys, imp

    try:
        __script__
    except NameError:
        raise ImportError('invalid invocation of psico.fullinit')

    imp.load_module('psico', None, os.path.dirname(__script__),
            ('', '', imp.PKG_DIRECTORY))

    import psico.fullinit

# vi: ts=4:sw=4:smarttab:expandtab
