

from importlib.util import find_spec

def _try_import(package_name):
    if find_spec(package_name):
        return True
    else:
        return False
