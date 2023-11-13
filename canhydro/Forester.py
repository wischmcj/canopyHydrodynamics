"""The workhorse class, leverages the others to get results"""

# import matplotlib.pyplot as plt
# from matplotlib.pyplot import cm
from __future__ import annotations

import os

from canhydro.CylinderCollection import CylinderCollection
from canhydro.global_vars import input_dir, log

NAME = "Forester"


# Class(es) intended to be the workhorse(s) that manages our objects


class CollectionManager:
    def __get__(self, obj, objtype):
        if obj is None:
            return Forester(objtype)
        else:
            raise AttributeError(f"Forester isn't accessible via {objtype} instances")


class Forester:
    #   Read in file names and create cylinder collection via CC class
    #
    #   Read in file, mostly an array of cylinders
    #
    #   Create graph
    #
    # initialize our object level variables for cylinder objects
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
        collections = []
        for file_obj in self.file_names:
            if file_name == None or file_obj.name == file_name:
                c = CylinderCollection()
                c.from_csv(file_obj, dir)
                collections.append(c)
        self.cylinder_collections = collections

    def networkSimplex():
        # can be used to calculate flows on graphs with demands
        # we could set a demand of generates X volume of flow
        print("nxs")
