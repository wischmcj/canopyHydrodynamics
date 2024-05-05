from types import ModuleType

def _try_import(name, src):
    try:
        name = None
        import_statment = f'from src import name'
        eval(import_statment)
        return True 
    except ImportError:
        return False