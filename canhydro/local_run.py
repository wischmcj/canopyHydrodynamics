import Cylinder
import utils

import global_vars as vars
import os 
from pathlib import Path
import shutil
import stat

DIR = vars.DIR

def on_rm_error( func, path, exc_info):
    # path contains the path of the file that couldn't be removed
    # let's just assume that it's read-only and unlink it.
    os.chmod( path, stat.S_IWRITE )
    os.unlink( path )

def create_dir_and_file(filename) -> None:
    print(type(filename))
    os.makedirs(filename,  exist_ok=True)
    f = open("test\demofile2.csv", "w")
    f.write("Now the file has more content!")
    f.close()

def del_dir(filename) -> None: 
    shutil.rmtree(filename, onerror = on_rm_error)

if __name__ == '__main__':
    file_path = ''.join([DIR,"test\\"])
    file = os.path.dirname(file_path)
    create_dir_and_file(file)
    print( Path(''.join([DIR,"test"])))
    utils.read_file_names( Path(''.join([DIR,"test"])))
    del_dir(file)
