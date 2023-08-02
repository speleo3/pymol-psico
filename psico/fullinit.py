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
    from pathlib import Path
    from importlib import util
    import sys

    try:
        _script_path = __script__  # noqa: F821 Undefined name
    except NameError:
        raise ImportError('invalid invocation of psico.fullinit') from None

    path = Path(_script_path).parent / "__init__.py"
    spec = util.spec_from_file_location('psico', path)
    module = util.module_from_spec(spec)
    sys.modules['psico'] = module
    spec.loader.exec_module(module)

    import psico.fullinit  # noqa: F401

# vi: ts=4:sw=4:smarttab:expandtab
