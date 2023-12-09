from __future__ import annotations

import csv
import os
import shutil
import stat
from typing import Union

import numpy as np
from numba import njit, prange

from canhydro.global_vars import input_dir, log, output_dir, time_stamp


@njit()
def stack(list_of_array, col: bool == True):
    """
    numba doesnt play well with np stacks, so I ha to do it myself
    """
    num_in = len(list_of_array)
    left_shape = list_of_array[0].shape[0]
    shape = (num_in, left_shape)
    stacked_array = np.empty(shape)
    for j in prange(len(list_of_array)):
        stacked_array[j] = list_of_array[j]
    return stacked_array if not col else stacked_array.T


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


def read_file_names(file_path=input_dir):
    """Reads in filenames to list"""
    paths = sorted(file_path.iterdir(), key=os.path.getmtime)
    fileNames = [f.name for f in paths if f.suffix == ".csv"]
    print(paths)
    print(fileNames)
    return fileNames


def save_file(
    file,
    out_file: Union(list[dict], dict),
    overwrite: bool = False,
    subdir: str = "agg",
    fileFormat=".csv",
    method="",
):
    dir = "/".join([str(output_dir), method, ""]).replace("/", "\\")
    ofname = "_".join([file, method]).replace("/", "\\")
    ofname_ext = "".join([ofname, fileFormat]).replace("/", "\\")
    aggname = "_".join(["agg", method, fileFormat]).replace("/", "\\")
    aggname_ext = "".join([aggname, fileFormat]).replace("/", "\\")
    folderExists = os.path.exists(dir)
    fileExists = os.path.exists(dir + ofname_ext)
    aggExists = os.path.exists(dir + aggname_ext)
    if not folderExists:
        os.makedirs(dir)

    to_write = []
    if isinstance(out_file, dict):
        out_file = [out_file]

    headers = list(out_file[0].keys())
    to_write.append(headers)
    log.info(f"{to_write}")
    for dic in out_file:
        cur_row = []
        for _, value in dic.items():
            cur_row.append(value)
        cur_row.append(time_stamp)
        to_write.append(cur_row)

    if fileExists:
        with open(dir + ofname_ext) as csv_file:
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
        log.info(f"{to_write}")
        with open(dir + ofname_ext, "w") as csv_file:
            writer = csv.writer(csv_file)
            for row in to_write:
                writer.writerow(row)


def intermitent_log(prog: int, whole: int, msg: str, freq: int = 0.0001):
    if np.random.uniform(0, 1, 1) < freq:
        log.info(msg + str(np.round((prog / whole) * 100, decimals=1)))
        print(msg + str(np.round((prog / whole) * 100, decimals=1)))


def lam_filter(objects, a_lambda: function, return_all: bool = False):
    """Takes in a lambda that filters on cylinder attrs"""
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


def count_lines(start, lines=0, header=True, begin_start=None):
    """
    Counts the lines contained in this project and prints results by
    file. For Vanity.
    """
    if header:
        print("{:>10} |{:>10} | {:<20}".format("ADDED", "TOTAL", "FILE"))
        print("{:->11}|{:->11}|{:->20}".format("", "", ""))

    for thing in os.listdir(start):
        thing = os.path.join(start, thing)
        if os.path.isfile(thing):
            if thing.endswith(".py"):
                with open(thing) as f:
                    newlines = f.readlines()
                    newlines = len(newlines)
                    lines += newlines

                    if begin_start is not None:
                        reldir_of_thing = "." + thing.replace(begin_start, "")
                    else:
                        reldir_of_thing = "." + thing.replace(start, "")

                    print(
                        "{:>10} |{:>10} | {:<20}".format(
                            newlines, lines, reldir_of_thing
                        )
                    )

    for thing in os.listdir(start):
        thing = os.path.join(start, thing)
        if os.path.isdir(thing):
            lines = countlines(thing, lines, header=False, begin_start=start)

    return lines


if __name__ == "__main__":
    count_lines(r"C:\Users\wisch\documents\gitprojects\canopyhydrodynamics")
