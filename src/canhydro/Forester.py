"""The workhorse class, leverages the others to get results"""

from __future__ import annotations

import os
import pickle

from src.canhydro.CylinderCollection import CylinderCollection
from src.canhydro.global_vars import input_dir, log, output_dir
from src.canhydro.utils import create_dir_and_file
NAME = "Forester"


# Class(es) intended to be the workhorse(s) that manages our objects


class CollectionManager:
    def __get__(self, obj, objtype):
        if obj is None:
            return Forester(objtype)
        else:
            raise AttributeError(f"Forester isn't accessible via {objtype} instances")


class Forester:
    def __init__(self, file_names="", directory=input_dir) -> None:
        self.file_names = file_names
        self.directory = directory
        self.cylinder_collections = []

    def get_file_names(self, dir=input_dir):
        #     os.chdir(''.join([vars.DIR,'input']))
        #     fullPath = Path(''.join([vars.DIR,'input']))
        log.info(f"Searching {dir} for files")
        paths = sorted(dir.iterdir(), key=os.path.getmtime)
        self.file_names = paths
        file_names = [f.name for f in paths if f.suffix == ".csv"]
        log.info(f"The following files found in {dir}: {file_names}")
        return paths

    def qsm_from_file_names(self, dir=input_dir, file_name: str = None):
        if self.file_names == "":
            self.get_file_names(dir)
        if file_name == None:
            log.error(
                "A file name must be provided. To scan all files in location, specify 'All'"
            )
            return
        collections = []
        for file_obj in self.file_names:
            if file_name == "All" or file_name in file_obj.name:
                c = CylinderCollection()
                c.from_csv(file_obj, dir)
                collections.append(c)
        if len(collections) == 0:
            log.error(f"File {file_name} not found in input directory {dir}")
        self.cylinder_collections = collections
