import __builtin__
from types import ModuleType


class DummyModule(ModuleType):
    def __getattr__(self, key):
        return None
    __all__ = []   # support wildcard imports


excluded_imports = ['geopandas']

def tryimport(name, globals={}, locals={}, fromlist=[], level=-1):
    try:
        if name in excluded_imports:
            return DummyModule(name)
        else:    
            return realimport(name, globals, locals, fromlist, level)
    except ImportError:
        return DummyModule(name)

realimport, __builtin__.__import__ = __builtin__.__import__, tryimport