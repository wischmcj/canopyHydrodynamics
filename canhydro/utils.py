from __future__ import annotations

import math
import os
import shutil
import stat
from typing import Union, Tuple, List
import warnings

from pandas import read_excel, ExcelWriter
import matplotlib.pyplot as plt 
from canhydro.global_vars import input_dir, log

# from global_vars import input_dir, log


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


def save_file(self, toWrite=[], subdir: str = "agg", fileFormat=".png", method=""):
    proj = "XY"
    if self.rev:
        proj = "XZ"
    file_arr = os.path.splitext(os.path.basename(self.filename))
    dir = "/".join([self.output_dir, method, ""]).replace("/", "\\")
    ofname = "_".join([file_arr[0], method, proj, fileFormat]).replace("/", "\\")
    aggname = "_".join(["agg", method, proj, fileFormat]).replace("/", "\\")
    folderExists = os.path.exists(dir)
    fileExists = os.path.exists(dir + ofname)
    aggExists = os.path.exists(dir + aggname)
    if not folderExists:
        os.makedirs(dir)
    if fileFormat == ".png":
        plt.savefig(dir + ofname, format="png", dpi=1200)
    else:
        if fileExists:
            exist = pd.read_excel(
                open(dir + ofname, "rb"), sheet_name=method, engine="openpyxl"
            )
            toWrite = toWrite.append(exist)
        with pd.ExcelWriter(dir + ofname, engine="openpyxl", mode="w") as writer:
            toWrite.to_excel(writer, index=False, sheet_name=method)
        if not aggExists:
            with pd.ExcelWriter(dir + aggname, engine="openpyxl", mode="w") as writer:
                toWrite.to_excel(writer, index=False, sheet_name=method)
        else:
            exist = pd.read_excel(
                open(dir + aggname, "rb"), sheet_name=method, engine="openpyxl"
            )
            toWrite = toWrite.append(exist)
            with pd.ExcelWriter(dir + aggname, engine="openpyxl", mode="w") as writer:
                toWrite.to_excel(writer, index=False, sheet_name=method)


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


def saveFile(self, toWrite=[], subdir: str = "agg", fileFormat=".png", method=""):
    proj = self._projection
    file_arr = os.path.splitext(os.path.basename(self._fileName))
    dir = "/".join([self._output_dir, method, ""]).replace("/", "\\")
    ofname = "_".join([file_arr[0], method, proj, fileFormat]).replace("/", "\\")
    aggname = "_".join(["agg", method, proj, fileFormat]).replace("/", "\\")
    folderExists = os.path.exists(dir)
    fileExists = os.path.exists(dir + ofname)
    aggExists = os.path.exists(dir + aggname)
    if not folderExists:
        os.makedirs(dir)
    if fileFormat == ".png":
        plt.savefig(dir + ofname, format="png", dpi=1200)
    else:
        if not fileExists:
            with pd.ExcelWriter(dir + ofname, engine="openpyxl", mode="w") as writer:
                toWrite.to_excel(writer, index=False, sheet_name=method)
        else:
            exist = pd.read_excel(
                open(dir + ofname, "rb"), sheet_name=method, engine="openpyxl"
            )
            toWrite = toWrite.append(exist)
            with pd.ExcelWriter(dir + ofname, engine="openpyxl", mode="w") as writer:
                toWrite.to_excel(writer, index=False, sheet_name=method)
        if not aggExists:
            with pd.ExcelWriter(dir + aggname, engine="openpyxl", mode="w") as writer:
                toWrite.to_excel(writer, index=False, sheet_name=method)
        else:
            exist = pd.read_excel(
                open(dir + aggname, "rb"), sheet_name=method, engine="openpyxl"
            )
            toWrite = toWrite.append(exist)
            with pd.ExcelWriter(dir + aggname, engine="openpyxl", mode="w") as writer:
                toWrite.to_excel(writer, index=False, sheet_name=method)

def countlines(start, lines=0, header=True, begin_start=None):
    """
    Counts the lines contained in this project and prints results by 
    file. For Vanity.
    """
    if header:
        print('{:>10} |{:>10} | {:<20}'.format('ADDED', 'TOTAL', 'FILE'))
        print('{:->11}|{:->11}|{:->20}'.format('', '', ''))

    for thing in os.listdir(start):
        thing = os.path.join(start, thing)
        if os.path.isfile(thing):
            if thing.endswith('.py'):
                with open(thing, 'r') as f:
                    newlines = f.readlines()
                    newlines = len(newlines)
                    lines += newlines

                    if begin_start is not None:
                        reldir_of_thing = '.' + thing.replace(begin_start, '')
                    else:
                        reldir_of_thing = '.' + thing.replace(start, '')

                    print('{:>10} |{:>10} | {:<20}'.format(
                            newlines, lines, reldir_of_thing))


    for thing in os.listdir(start):
        thing = os.path.join(start, thing)
        if os.path.isdir(thing):
            lines = countlines(thing, lines, header=False, begin_start=start)

    return lines

# if __name__ == "__main__":
#     countlines(r'C:\Users\wisch\documents\gitprojects\canopyhydrodynamics')