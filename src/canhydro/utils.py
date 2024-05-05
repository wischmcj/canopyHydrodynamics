from __future__ import annotations

import csv
import os
import shutil
import stat

import numpy as np
from src.canhydro.global_vars import input_dir, log, output_dir, time_stamp

def on_rm_error(func, path, exc_info):
    # path contains the path of the file that couldn't be removed
    # let's just assume that it's read-only and unlink it.
    os.chmod(path, stat.S_IWRITE)
    os.unlink(path)


def create_dir_and_file(filename,) -> None:
    print(type(filename))
    os.makedirs(filename, exist_ok=True)
    # f = open(filename, "w+")
    # f.write("Now the file has more content!")
    # f.close()


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
    out_file,
    overwrite: bool = False,
    subdir: str = "agg",
    fileFormat=".csv",
    method=""
):
    dir = "/".join([str(output_dir), method, ""])
    ofname = "_".join([file, method])
    ofname_ext = "".join([ofname, fileFormat])
    aggname = "_".join(["agg", method, fileFormat])
    aggname_ext = "".join([aggname, fileFormat])
    folderExists = os.path.exists(dir)
    fileExists = os.path.exists(dir + ofname_ext)
    aggExists = os.path.exists(dir + aggname_ext)
    if not folderExists:
        log.info(f'folder doesnt exist, creating {dir}')
        os.makedirs(dir)
    to_write = []
    if isinstance(out_file, dict):
        out_file = [out_file]

    headers = list(out_file[0].keys())
    log.info(f'file headers {headers}')
    to_write.append(headers)
    for dic in out_file:
        cur_row = []
        for _, value in dic.items():
            cur_row.append(value)
        cur_row.append(time_stamp)
        to_write.append(cur_row)
    # breakpoint()
    # if fileExists:
    #     with open(dir + ofname_ext) as csv_file:
    #         reader = csv.reader(csv_file)
    #         existing_rows = [x for x in reader]
    #         if existing_rows[0] == headers:
    #             for row in existing_rows[1:]:
    #                 if row != []:
    #                     to_write.append(row)
    #         else:
    #             log.warning(
    #                 f"Existing {ofname_ext} file has different headers, to overwrite pass overwrite =true"
    #             )
    if overwrite:
        # log.info(f"{to_write}")
        file = dir + ofname_ext
        log.info(f"attempting to write to {file}")
        with open(file, "w") as csv_file:
            writer = csv.writer(csv_file)
            for row in to_write:
                writer.writerow(row)
    return file


def intermitent_log(prog: int, whole: int, msg: str, freq: int = 0.0001):
    if np.random.uniform(0, 1, 1) < freq:
        log.info(msg + str(np.round((prog / whole) * 100, decimals=1)))
        print(msg + str(np.round((prog / whole) * 100, decimals=1)))


def lam_filter(objects, a_lambda: function, return_all: bool = False):
    """Takes in a lambda that filters on cylinder attrs
        returns a list of cylinders meeting that requirements 
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
