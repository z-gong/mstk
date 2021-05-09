_OPENBABEL_IMPORTED = False
try:
    import pybel
except ImportError:
    pass
else:
    import openbabel
    _OPENBABEL_IMPORTED = True

if not _OPENBABEL_IMPORTED:
    try:
        from openbabel import pybel
    except ImportError:
        pass
    else:
        from openbabel import openbabel

        _OPENBABEL_IMPORTED = True

if not _OPENBABEL_IMPORTED:
    raise Exception('OpenBabel not found')
