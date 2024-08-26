from __future__ import annotations

import logging
import os
from pathlib import Path

from canopyhydro.configuration import input_dir
from canopyhydro.CylinderCollection import CylinderCollection

log = logging.getLogger("model")


NAME = "Forester"


class Forester:
    def __init__(self, directory=input_dir) -> None:
        self.file_names = []
        self.directory = directory
        self.cylinder_collections = []
        self.file_names = self.get_file_names()

    def get_file_names(self):
        directory = self.directory
        log.debug(f"Searching {directory} for files")
        directory = Path(directory)
        paths = sorted(directory.iterdir(), key=os.path.getmtime)
        self.file_names = paths
        file_names = [f.name for f in paths if f.suffix == ".csv"]
        log.debug(f"The following files found in {directory}: {file_names}")
        return paths

    def qsm_to_collection(self, file_name: str = "All", directory: Path = ""):
        """Creates a Cylinder collection from the given file
            and adds the collection to the Forester object

        Args:
            file_name (str):
                Name of the csv file containing QSM data for
                the desired cylinder collection
            directory (str, optional):
                File directory in which the specified file is located.
                Defaults to input_dir.
        """
        if directory == "":
            directory = self.directory
        if self.file_names == "":
            self.get_file_names()
        collections = []
        for file_obj in self.file_names:
            if ".csv" not in file_name and file_name != "All":
                file_name = file_name + ".csv"
            if file_name == "All" or file_obj.name == file_name:
                c = CylinderCollection()
                c.from_csv(file_obj)
                collections.append(c)
        if len(collections) == 0:
            msg = f"File {file_name} not found in input directory {directory}"
            log.error(msg)
            raise FileNotFoundError(msg)
        self.cylinder_collections = collections
