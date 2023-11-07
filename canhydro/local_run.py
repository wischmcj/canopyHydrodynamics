from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path

import global_vars
from Forester import Forester
import utils

DIR = global_vars.DIR
test_input_dir = global_vars.test_input_dir


def on_rm_error(func, path, exc_info):
    # path contains the path of the file that couldn't be removed
    # let's just assume that it's read-only and unlink it.
    os.chmod(path, stat.S_IWRITE)
    os.unlink(path)


def create_dir_and_file(filename) -> None:
    print(type(filename))
    os.makedirs(filename, exist_ok=True)
    f = open(r"test\demofile2.csv", "w")
    f.write("Now the file has more content!")
    f.close()


def del_dir(filename) -> None:
    shutil.rmtree(filename, onerror=on_rm_error)


if __name__ == "__main__":
    forest = Forester()
    forest.qsm_from_file_names(test_input_dir)
