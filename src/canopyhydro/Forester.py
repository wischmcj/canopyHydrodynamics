from __future__ import annotations

import logging
import os
from pathlib import Path

import toml
from src.canhydro.CylinderCollection import CylinderCollection

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    input_dir = Path(config["directories"]["input_dir"])

log = logging.getLogger("model")


NAME = "Forester"


class Forester:
    def __init__(self, file_names="", directory=input_dir) -> None:
        self.file_names = file_names
        self.directory = directory
        self.cylinder_collections = []

    def get_file_names(self, dir=input_dir):
        log.info(f"Searching {dir} for files")
        dir = Path(dir)
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
            if ".csv" not in file_name:
                file_name = file_name + ".csv"
            if file_name == "All" or file_obj.name == file_name:
                c = CylinderCollection()
                c.from_csv(file_obj, dir)
                collections.append(c)
        if len(collections) == 0:
            msg = f"File {file_name} not found in input directory {dir}"
            log.error(msg)
            raise FileNotFoundError(msg)
        self.cylinder_collections = collections
