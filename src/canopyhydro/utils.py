from __future__ import annotations

import calendar
import csv
import logging
import os
import shutil
import stat
import time
from importlib.util import find_spec
from typing import Callable

import numpy as np

from canopyhydro.configuration import input_dir, output_dir

log = logging.getLogger("model")

current_GMT = time.gmtime()
time_stamp = str(calendar.timegm(current_GMT))


def _try_import(package_name):
    if find_spec(package_name):
        return True
    else:
        return False


has_numba = _try_import("numba")
if has_numba:
    from numba import njit, prange
    from numba.typed import List

# Data munging utils


def stack(to_stack: list[np.array], col: bool = True):
    """
    A wrapper for njit stack that handles errors and allows for
    less strict typing
    """
    if not has_numba:
        ret = non_njit_stack(to_stack, col)
    else:
        list_of_array = List(to_stack)
        try:
            ret = njit_stack(list_of_array, col)
        except ValueError:
            ret = non_njit_stack(to_stack, col)
    return ret


def non_njit_stack(to_stack: list[np.array], col: bool = True):
    list_of_array = list(to_stack)
    num_in = len(list_of_array)
    left_shape = list_of_array[0].shape[0]
    shape = (num_in, left_shape)
    stacked_array = np.empty(shape)
    for idx, _ in enumerate(list_of_array):
        stacked_array[idx] = list_of_array[idx]
    return stacked_array if not col else stacked_array.T


if has_numba:

    @njit()
    def njit_stack(list_of_array: np.array[np.array()], col: bool):
        """
        numba doesn't play well with np stacks, so I had to do it myself
        """
        num_in = len(list_of_array)
        left_shape = list_of_array[0].shape[0]
        shape = (num_in, left_shape)
        stacked_array = np.empty(shape)
        for j in prange(len(list_of_array)):
            stacked_array[j] = list_of_array[j]
        return stacked_array if not col else stacked_array.T
# File system utils


def on_rm_error(path):
    os.chmod(path, stat.S_IWRITE)
    os.unlink(path)


def create_dir_and_file(
    filename,
) -> None:
    os.makedirs(filename, exist_ok=True)


def del_dir(filename) -> None:
    shutil.rmtree(filename, onerror=on_rm_error)


def read_file_names(file_path=input_dir):
    """Reads in filenames to list"""
    paths = sorted(file_path.iterdir(), key=os.path.getmtime)
    fileNames = [f.name for f in paths if f.suffix == ".csv"]
    log.info(f"fileNames found in {file_path} on read: {fileNames}")
    return fileNames


def save_file(
    file: str,
    out_file: dict | list[dict],
    overwrite: bool = False,
    fileFormat: str = ".csv",
    method: str = "",
):
    """
    A somewhat overly complex file saving function.
    If file exists and overwrite = false, it will append to the file

    Args:
        file: str
            The name of the file to write to
        out_file: Union[dict, list[dict]]
            The data to write to the file
        overwrite: bool
            Whether or not to overwrite the file if it exists
        fileFormat: str
            The file format to save the file as ()
        method: str
    Note:
        'agg*' variables intended to support future
            updates - adding additional write to an
            append only file for better history tracking
    """
    dir = "/".join([str(output_dir), method, ""])
    ofname = "_".join([file, method])
    ofname_ext = "".join([ofname, fileFormat])
    folderExists = os.path.exists(dir)
    fileExists = os.path.exists(dir + ofname_ext)

    # aggname = "_".join(["agg", method, fileFormat])
    # aggname_ext = "".join([aggname, fileFormat])
    # aggExists = os.path.exists(dir + aggname_ext)#

    if not folderExists:
        os.makedirs(dir)

    to_write = []
    if isinstance(out_file, dict):
        out_file = [out_file]

    headers = list(out_file[0].keys())
    to_write.append(headers)
    log.info(f"{to_write}")
    # adding each dict in out_file to a
    #   single to_write list[list]
    for dic in out_file:
        cur_row = []
        for _, value in dic.items():
            cur_row.append(value)
        cur_row.append(time_stamp)
        to_write.append(cur_row)

    if fileExists and not overwrite:
        # Reading existing data into to_write in
        #  order to append new data
        with open(dir + ofname_ext, "w+") as csv_file:
            reader = csv.reader(csv_file)
            existing_rows = list(reader)
            if existing_rows[0] == headers:
                for row in existing_rows[1:]:
                    if row != []:
                        to_write.append(row)
            else:
                log.warning(
                    f"Existing { ofname_ext} file has different headers, to overwrite pass ovewrite =true"
                )
    if overwrite:
        # log.info(f"{to_write}")
        file = dir + ofname_ext
        log.info(f"attempting to write to {file}")
        with open(file, "w") as csv_file:
            writer = csv.writer(csv_file)
            for row in to_write:
                writer.writerow(row)
    return file


# other utils


def intermitent_log(prog: int, whole: int, msg: str, freq: int = 0.0001):
    if np.random.uniform(0, 1, 1) < freq:
        log.info(msg + str(np.round((prog / whole) * 100, decimals=1)))
        print(msg + str(np.round((prog / whole) * 100, decimals=1)))


def lam_filter(objects, a_lambda: Callable, return_all: bool = False):
    """
    Takes in a lambda that filters on cylinder attrs
    returns a list of cylinders for which that
    lambda func returns true
    
    Args:
        return_all (bool): if true, returns all cylinders along with
            a boolean array of the same length w/ the
            results of the lambda function passed

    Example:
        lam_filter(cylinders, lambda: diameter > 0.5)
        returns: all cylinders with diameter > 0.5
    """
    if a_lambda.__code__.co_name != "<lambda>":
        raise ValueError("a lambda required")
    if return_all:
        filtered = [(obj, eval(a_lambda.__code__, vars(obj).copy())) for obj in objects]
    else:
        filtered = [
            (obj, True) for obj in objects if eval(a_lambda.__code__, vars(obj).copy())
        ]
    objs, res = zip(*filtered)
    return objs, res
