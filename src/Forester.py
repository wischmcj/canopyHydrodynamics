"""The workhorse class, leverages the others to get results"""

# import matplotlib.pyplot as plt
# from matplotlib.pyplot import cm
from pathlib import Path
# from random import random
# from multiprocessing import Pool

# from shapely.geometry import Polygon, Point
# from shapely.ops import unary_union, transform
# from descartes import PolygonPatch
# from mpl_toolkits import mplot3d
# from pickle import dump, load


from time import sleep

# import networkx as nx
# import openpyxl
# import geopandas as geo
# import numpy as np
import calendar
# import time
# import copy
# import math
import logging
# import settings
import os

import global_vars as vars


from pandas import to_excel as pd

time_stamp = str(calendar.timegm(current_GMT))


NAME = "Cylinder"
# logging.basicConfig(filename=''.join(['log_',str(time_stamp)])  , filemode='w', level=logging.DEBUG, encoding='utf-8',level=os.environ.get("LOGLEVEL", "INFO"))
log = logging.getLogger("my-logger")
 
#Class intented to be the workhorse that manages our objects
class Forester:
    #   Read in file names and create cylinder collection via CC class
    # 
    #   Read in file, mostly an array of cylinders 
    #   
    #   Create graph
    #       

    #initialize our object level variables for cylider objects 
    def __init__(self, file_names = np.nan, directory = DIR) -> None:
        self._file_names = file_names
        self._directory =directory
        self.__cylinder_collections = []
    
    def getFileNames(self):
    #     os.chdir(''.join([vars.DIR,'input']))
    #     fullPath = Path(''.join([vars.DIR,'input']))
        os.chdir(''.join([self._directory,'input']))
        full_path = Path(''.join([self._directory,'input']))
        paths = sorted(full_path.iterdir(),key=os.path.getmtime)
        file_names = [f.name for f in paths if  f.suffix == '.csv' ]
        self._file_names = file_names
        return file_names
    
    def QSM_FromFileNames(self):
        collections = []
        for filename in self._file_names:
            c = CylinderCollection(filename, directory = DIR)
            c.read_csv()
            collections.append(c)

        return collections

    def networkSimplex():
        #can be used to calculate flows on graphs with demands 
        #we could set a demand of generates X volume of flow 
        print('nxs')