from __future__ import annotations

import os
import stat


def within_range(expected, actual, err):
    return actual >= expected - expected * err and actual <= expected + expected * err


def on_rm_error(func, path, exc_info):
    # path contains the path of the file that couldn't be removed
    # let's just assume that it's read-only and unlink it.
    os.chmod(path, stat.S_IWRITE)
    os.unlink(path)


def create_dir_and_file(filename) -> None:
    print(type(filename))
    os.makedirs(filename, exist_ok=True)
    f = open("test\\demofile2.csv", "w")
    f.write("Now the file has more content!")
    f.close()
