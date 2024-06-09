
import os
import toml
import sys

sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())

import logging

log = logging.getLogger("model")

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    test_input_dir = config["directories"]['test_input_dir']

from src.canhydro.CylinderCollection import unpickle_collection, pickle_collection
from src.canhydro.Forester import Forester



def initialize_collection(file = "5_SmallTree", from_pickle = False, **kwargs):
    if from_pickle:
        collection = load_from_pickle(**kwargs)
    else:
        log.info(f"initializing collection...")
        try:
            forest = Forester()
            forest.get_file_names(dir=test_input_dir)
            forest.qsm_from_file_names(file_name=file)
            collection = forest.cylinder_collections[0]
            collection.project_cylinders("XY")
        except Exception as e:
            log.info(f"Error initializing collection for file {file}: {e}")
            return None
        log.info(f"successfully initialized collection")
        # try_pickle_collection(basic_collection,file)
    return collection

def load_from_pickle(file, pickle_point, angle):
    pickle_file = f'{file}_pickle__{pickle_point}_{angle}'
    collection = None
    try:
        log.info(f'Loading collection for {file}, case {angle} from pickle file {pickle_file}')
        collection = unpickle_collection(pickle_file)  
        if collection: 
            log.info(f'Collection for {file}, case {angle} successfully loaded from pickle file {pickle_file}')
    except Exception as e:
        log.info(f'failed to load pickle for {file}, case {angle}  :{e}')
    return collection

def try_pickle_collection(collection, designation = ""):
    log.info(f"attempting to pickle for file {collection.file_name}, des {designation}")
    # pickle_file = pickle_collection(collection, designation)
    try: 
        pickle_file = pickle_collection(collection, designation)
        log.info(f"successfully created pickle for file {collection.file_name}")
    except Exception as e:
        log.info(f"Error pickling file {collection.file_name}: {e}")
        return
