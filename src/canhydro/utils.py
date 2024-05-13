from __future__ import annotations

import csv
import calendar
import os
import shutil
import stat
import logging
import toml
import time 
from typing import Union
from pathlib import Path

import numpy as np

from src.canhydro.import_options import _try_import

has_numba = _try_import('numba')
if has_numba:
    from numba import njit, prange
    from numba.typed import List

log = logging.getLogger("model")

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    input_dir = Path(config["directories"]['input_dir'])
    output_dir = Path(config["directories"]['output_dir'])
    
 
current_GMT = time.gmtime()
time_stamp = str(calendar.timegm(current_GMT))


# Data munging utils

def stack(to_stack:list[np.array], col: bool = True):
    """
        A wrapper for njit stack that handles errors and allows for 
        less strict typing 
    """
    if not has_numba:
        non_njit_stack(to_stack, col)
    else:
        list_of_array = List(to_stack)
        try:
            njit_stack(list_of_array,col)
        except ValueError as err:
            non_njit_stack(to_stack, col)


def non_njit_stack(to_stack:list[np.array], col: bool = True):
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
    def njit_stack(list_of_array:np.array[np.array()], col: bool):
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


def create_dir_and_file(filename,) -> None:
    print(type(filename))
    os.makedirs(filename, exist_ok=True)


def del_dir(filename) -> None:
    shutil.rmtree(filename, onerror=on_rm_error)


def read_file_names(file_path=input_dir):
    """Reads in filenames to list"""
    paths = sorted(file_path.iterdir(), key=os.path.getmtime)
    fileNames = [f.name for f in paths if f.suffix == ".csv"]
    log.info(f'fileNames found in {file_path} on read: {fileNames}')
    return fileNames


def save_file(
    file: str,
    out_file: Union[dict, list[dict]],
    overwrite: bool = False,
    fileFormat: str = ".csv",
    method: str ="",
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
<<<<<<< HEAD
=======

    # aggname = "_".join(["agg", method, fileFormat])
    # aggname_ext = "".join([aggname, fileFormat])
    # aggExists = os.path.exists(dir + aggname_ext)#

    if not folderExists:
        os.makedirs(dir)
>>>>>>> 8aeaf6b (Revert "Revert "Merge branch 'main' into ingest-validation-data"")

    # aggname = "_".join(["agg", method, fileFormat])
    # aggname_ext = "".join([aggname, fileFormat])
    # aggExists = os.path.exists(dir + aggname_ext)#

    if not folderExists:
        log.info(f'folder doesnt exist, creating {dir}')
        os.makedirs(dir)
    to_write = []
    if isinstance(out_file, dict):
        out_file = [out_file]

    headers = list(out_file[0].keys())
    log.info(f'file headers {headers}')
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
    log.info(f"{to_write}")
    with open(dir + ofname_ext, "w") as csv_file:
        writer = csv.writer(csv_file)
        for row in to_write:
            writer.writerow(row)

# other utils

def intermitent_log(prog: int, whole: int, msg: str, freq: int = 0.0001):
    if np.random.uniform(0, 1, 1) < freq:
        log.info(msg + str(np.round((prog / whole) * 100, decimals=1)))
        print(msg + str(np.round((prog / whole) * 100, decimals=1)))


def lam_filter(objects, a_lambda: function, return_all: bool = False):
    """
        Takes in a lambda that filters on cylinder attrs
        returns a list of cylinders for which that
        lambda func returns true
        i.e. lam_filter(cylinders, lambda: diameter > 0.5)
                -returns all cylinders with diameter > 0.5

        return_all: bool 
            - if true, returns all cylinders along with 
                a boolean array of the same length w/ the 
                results of the lambda function passed
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
