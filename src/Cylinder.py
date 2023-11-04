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

import global_vars


from pandas import to_excel as pd

time_stamp = str(calendar.timegm(current_GMT))

NAME = "Cylinder"
logging.basicConfig(filename=''.join(['log_',str(time_stamp)])  , filemode='w', level=logging.DEBUG, encoding='utf-8',level=os.environ.get("LOGLEVEL", "INFO"))
log = logging.getLogger("my-logger")
 
class CylinderCollection:
    # A cylinder collection is just an array of cylinder objects 
    #  with some of its own variables 
    # Its purpose is to act as a storage for entire trees as well as 
    # sub trees. This might include the stem flow, drip flow
    # 
    #  

    #initialize our object level variables for cylider objects 
    def __init__(self, CylinderCollection) -> None:
        #Base Attributes from file read in
        # self.cylinderCollection = CylinderCollection
        self.dimensions = {
                            "x"                : np.nan,
                            "y"                : np.nan,
                            "z"                : np.nan,
                            "radius"           : np.nan,
                            "length"           : np.nan,
                            }
        self.model_loc = {
                            "branch_order"     : np.nan,
                            "branch_id"        : np.nan,
                            "volume"           : np.nan,
                            "parent_id"        : np.nan,
                            "rev_branch_order" : np.nan,
                            "section_id"       : np.nan
                        }
        self.x                  = np.nan #len 2 array
        self.y                  = np.nan #len 2 array
        self.z                  = np.nan #len 2 array
        self.radius             = np.nan
        self.length             = np.nan
        self.branch_order       = np.nan
        self.branch_id          = np.nan
        self.volume             = np.nan
        self.parent_id          = np.nan
        self.rev_branch_order   = np.nan
        self.section_id         = np.nan
        
        #Calculated based off of projection
        self.projected_data= {
                'XY':{'unit_vect':[],
                        'rev_unit_vect':[],
                        'rise_angle_rad':-7,
                        'polygon': None},
                'YZ':{'unit_vect':[],
                        'rev_unit_vect':[],
                        'rise_angle_rad':-7,
                        'polygon': None}
        }

        #graph attributes, all in 3-D
        self.graph                  = nx.Graph()
        self.flow_id                = np.nan
        self.flow_type              = np.nan
        self.begins_at_drip_point   = np.nan #bool
        self.begins_at_divide_point = np.nan #bool

        #others
        self.stem_path_id = np.nan
        
    def get_flow_data():
        #Replace trunk with single node 

        #Remove edges with out flow 
        #Find the connected component with root in it, this is stem flow 

        #Remove all stem flow from orig graph, these are drip paths
            #Find node of minimal height in each component, this is the drip point 
    
    def 