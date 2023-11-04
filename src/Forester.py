"""The workhorse class, leverages the others to get results"""

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from pathlib import Path
from random import random
from multiprocessing import Pool

from shapely.geometry import Polygon, Point
from shapely.ops import unary_union, transform
from descartes import PolygonPatch
from mpl_toolkits import mplot3d
from pickle import dump, load


from time import sleep

import networkx as nx
import openpyxl
import geopandas as geo
import numpy as np
import calendar
import time
import copy
import math
import logging
import settings
import os

import global_vars as vars


from pandas import to_excel as pd

time_stamp = str(calendar.timegm(current_GMT))

NAME = "Cylinder"
logging.basicConfig(filename=''.join(['log_',str(time_stamp)])  , filemode='w', level=logging.DEBUG, encoding='utf-8',level=os.environ.get("LOGLEVEL", "INFO"))
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
    def __init__(self, filename) -> None:
        self.variable = v
    
    def get_file_names():
        os.chdir(''.join([vars.DIR,'input']))
        fullPath = Path(''.join([vars.DIR,'input']))
        paths = sorted(fullPath.iterdir(),key=os.path.getmtime)
        fileNames = [f.name for f in paths if  f.suffix == '.csv' ]
        print(fileNames)
        return fileNames
   
    def network_simplex():
        #can be used to calculate flows on graphs with demands 
        #we could set a demand of generates X volume of flow 