"""These constitute our way of 'copying' the QSMs"""


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

    #initialize our object level variables for cylider objects 
    def __init__(self, filename) -> None:
        self.filename               = filename
        self.contained_cylinders              = np.nan #np array of cyls
        # self.df = pd.DataFrame() # Full read in data frame?
        
        # self.noCylinders = np.nan #Just the len of contained_cylinders
        
        #Aggregate values from file
        self.component_surface_area              = np.nan 
        self.component_volume              = np.nan
        self.max_branch_order              = np.nan
        self.max_rev_branch_order             = np.nan
        self.canopy_scope                          = np.nan #desc of canopy
        self.extent                           = {'min':[np.nan,np.nan,np.nan],'max':[np.nan,np.nan,np.nan]}
            # to populate with x,y,z mins and maxs
        self.aggregate_angle              = np.nan
        self.descriptive_vectors              = np.nan #Average, median, mode vectors
        self.treeQualities               = pd.DataFrame({'total_psa':-1 ,
                                            'tot_hull_area':-1, 
                                            'stem_flow_hull_area':-1 ,
                                            'stem_psa':-1 ,
                                            'flowStats':-1 ,
                                            'DBH':-1,
                                            'tot_surface_area':-1,
                                            'stem_surface_area':-1
                                            },index=[0])

        #Projection Attrs
        self.union_poly              = None
        self.stem_path_lengths              = []
        self.hull                          = np.nan
        self.stem_hull             = np.nan

        #Special case tree attributes
        self.stem_paths              =  [[]] #Cyl collection?
        self.trunk                           = [] #Collection of cylinders? id list?
        
        #Graph and Attributes        
        self.graph                           = nx.Graph()

        self.flows                          =[{
                                                 'cyls':[],
                                                 'drip_point':id,
                                                 'attributes':
                                                         {'cyls' : 0, 'len' : 0, 'sa' : 0, 'pa' : 0, 'as' : 0},
                                             }]
        self.drip_points              = {'x':np.nan, 'y':np.nan, 'z':np.nan, 'flow_id':np.nan}
        self.flow_to_drip             = {0:1} # A dictionary of flow ids with values equal to their drip node ids 
        self.trunk_nodes              = []
        self.drip_loc              = np.nan

        #Calculations using graph results
        self.stemTotal              = {
                                         'attributes':
                                             {'cyls' : 0, 'len' : 0, 'sa' : 0, 'pa' : 0, 'as' : 0},
                                         'loc': {'x':np.nan, 'y':np.nan, 'z':np.nan}
                                     }
        self.divide_points              = []
        self.stemPolys              = []
        self.compGraphs              = []

    def save_file(self, toWrite = [], subdir:str = 'agg', fileFormat ='.png',method=''):
        proj = 'XY'
        if self.rev : 
            proj='XZ'
        file_arr = os.path.splitext(os.path.basename(self.filename))
        dir = '/'.join([self.output_dir, method, '']).replace('/','\\')
        ofname = '_'.join([file_arr[0], method, proj , fileFormat ]).replace('/','\\')
        aggname = '_'.join(['agg', method, proj , fileFormat ]).replace('/','\\')
        folderExists = os.path.exists(dir)
        fileExists = os.path.exists(dir+ofname)
        aggExists = os.path.exists(dir+aggname)
        if not folderExists:
            os.makedirs(dir)
        if fileFormat =='.png': 
            plt.savefig(dir+ofname, format='png', dpi=1200)
        else:
            if fileExists:
                exist = pd.read_excel(open(dir+ofname, 'rb'), sheet_name=method, engine = 'openpyxl')  
                toWrite = toWrite.append(exist)
            with pd.ExcelWriter(dir + ofname, engine = 'openpyxl',mode='w') as writer:
                toWrite.to_excel(writer, index = False, sheet_name =method)
            if not aggExists:
                with pd.ExcelWriter(dir + aggname, engine = 'openpyxl',mode='w') as writer:
                    toWrite.to_excel(writer, index = False, sheet_name =method)    
            else:
                exist = pd.read_excel(open(dir+aggname, 'rb'), sheet_name=method, engine = 'openpyxl')  
                toWrite = toWrite.append(exist)
                with pd.ExcelWriter(dir + aggname, engine = 'openpyxl',mode='w') as writer:
                    toWrite.to_excel(writer, index = False, sheet_name =method)  
    def read_csv(self,df = pd.DataFrame(), polys = [],projection='XY'):
        #Columns [ID?,ParentID?,x1,y1,z1,x2,y2,z2,radius,?,?,lenght,? ,? ,? ,? ,? ,? ,? ,BO]
        #Colnums [1  ,2        ,3 , 4, 5, 6, 7, 8,9    ,10,11,12   ,13,14,15,16,17,18,19,20]
        # x = x2-x1, y =y2-y1, z=z2-z1
        #number of cyliders = cnt values for radius
        # Theta  = angle made buy the cylinder axis
        self.projection =projection 
        if('partial' not in self.filename):
            self.df_full = pd.read_csv(self.filename,header=0)
        else: 
            self.df_full = df
        self.maxBO = np.max(self.BO)
        #self.df = self.df_full 
        self.df=self.df_full#.iloc[:130,:]
        #columns 3 and 6 represent our x values
        if projection == 'XZ':
            self.x = np.transpose(self.df.iloc[:,[3,6]].to_numpy())
            self.y = np.transpose(self.df.iloc[:,[5,8]].to_numpy())
            self.z = np.transpose(self.df.iloc[:,[4,7]].to_numpy())
        elif projection == 'YZ':
            self.x = np.transpose(self.df.iloc[:,[5,8]].to_numpy())
            self.y = np.transpose(self.df.iloc[:,[3,6]].to_numpy())
            self.z = np.transpose(self.df.iloc[:,[4,7]].to_numpy())     
        else: # 'XY'
            self.x = np.transpose(self.df.iloc[:,[3,6]].to_numpy())   
            self.y = np.transpose(self.df.iloc[:,[4,7]].to_numpy()) 
            self.z = np.transpose(self.df.iloc[:,[5,8]].to_numpy()) 
        # for side view
        self.cylID = self.df.iloc[:,1].to_numpy()
        self.pID =self.df.iloc[:,2].to_numpy()
        self.radius = self.df.iloc[:,9].to_numpy()
        self.noCylinders = self.radius.size
        self.cLength = self.df.iloc[:,12].to_numpy()
        self.BO = self.df.iloc[:,20].to_numpy()
        self.branchID = self.df.iloc[:,4].to_numpy()
        self.maxBO = np.max(self.BO)
        self.bID = self.df.iloc[:,24].to_numpy()
        if projection == 'XZ':
            self.dx = self.df.iloc[:,6].to_numpy() - self.df.iloc[:,3].to_numpy()
            self.dy = self.df.iloc[:,8].to_numpy() - self.df.iloc[:,5].to_numpy()
            self.dz = self.df.iloc[:,7].to_numpy() - self.df.iloc[:,4].to_numpy()
        elif projection == 'YZ':
            self.dx = self.df.iloc[:,8].to_numpy() - self.df.iloc[:,5].to_numpy()
            self.dy = self.df.iloc[:,6].to_numpy() - self.df.iloc[:,3].to_numpy()
            self.dz = self.df.iloc[:,7].to_numpy() - self.df.iloc[:,4].to_numpy()    
        else: # 'XY'
            self.dx = self.df.iloc[:,6].to_numpy() - self.df.iloc[:,3].to_numpy()  
            self.dy = self.df.iloc[:,7].to_numpy() - self.df.iloc[:,4].to_numpy()
            self.dz = self.df.iloc[:,8].to_numpy() - self.df.iloc[:,5].to_numpy() 
        if('partial' in self.filename): self.pSV = polys
        self.theta = np.arctan(self.dz/np.sqrt(self.dx**2 + self.dy**2))
        self.output_dir = ''.join([DIR,'output/'])
        log.info(self.filename + " initialized")

    def read_csv(self,df = pd.DataFrame(), polys = [],projection='XY'):
        #Columns [ID?,ParentID?,x1,y1,z1,x2,y2,z2,radius,?,?,lenght,? ,? ,? ,? ,? ,? ,? ,BO]
        #Colnums [1  ,2        ,3 , 4, 5, 6, 7, 8,9    ,10,11,12   ,13,14,15,16,17,18,19,20]
        # x = x2-x1, y =y2-y1, z=z2-z1
        #number of cyliders = cnt values for radius
        # Theta  = angle made buy the cylinder axis
        self.projection =projection 
        if('partial' not in self.filename):
            self.df_full = pd.read_csv(self.filename,header=0)
        else: 
            self.df_full = df
        self.maxBO = np.max(self.BO)
        #self.df = self.df_full 
        self.df=self.df_full#.iloc[:130,:]
        #columns 3 and 6 represent our x values
        if projection == 'XZ':
            self.x = np.transpose(self.df.iloc[:,[3,6]].to_numpy())
            self.y = np.transpose(self.df.iloc[:,[5,8]].to_numpy())
            self.z = np.transpose(self.df.iloc[:,[4,7]].to_numpy())
        elif projection == 'YZ':
            self.x = np.transpose(self.df.iloc[:,[5,8]].to_numpy())
            self.y = np.transpose(self.df.iloc[:,[3,6]].to_numpy())
            self.z = np.transpose(self.df.iloc[:,[4,7]].to_numpy())     
        else: # 'XY'
            self.x = np.transpose(self.df.iloc[:,[3,6]].to_numpy())   
            self.y = np.transpose(self.df.iloc[:,[4,7]].to_numpy()) 
            self.z = np.transpose(self.df.iloc[:,[5,8]].to_numpy()) 
        # for side view
        self.cylID = self.df.iloc[:,1].to_numpy()
        self.pID =self.df.iloc[:,2].to_numpy()
        self.radius = self.df.iloc[:,9].to_numpy()
        self.noCylinders = self.radius.size
        self.cLength = self.df.iloc[:,12].to_numpy()
        self.BO = self.df.iloc[:,20].to_numpy()
        self.branchID = self.df.iloc[:,4].to_numpy()
        self.maxBO = np.max(self.BO)
        self.bID = self.df.iloc[:,24].to_numpy()
        if projection == 'XZ':
            self.dx = self.df.iloc[:,6].to_numpy() - self.df.iloc[:,3].to_numpy()
            self.dy = self.df.iloc[:,8].to_numpy() - self.df.iloc[:,5].to_numpy()
            self.dz = self.df.iloc[:,7].to_numpy() - self.df.iloc[:,4].to_numpy()
        elif projection == 'YZ':
            self.dx = self.df.iloc[:,8].to_numpy() - self.df.iloc[:,5].to_numpy()
            self.dy = self.df.iloc[:,6].to_numpy() - self.df.iloc[:,3].to_numpy()
            self.dz = self.df.iloc[:,7].to_numpy() - self.df.iloc[:,4].to_numpy()    
        else: # 'XY'
            self.dx = self.df.iloc[:,6].to_numpy() - self.df.iloc[:,3].to_numpy()  
            self.dy = self.df.iloc[:,7].to_numpy() - self.df.iloc[:,4].to_numpy()
            self.dz = self.df.iloc[:,8].to_numpy() - self.df.iloc[:,5].to_numpy() 
        if('partial' in self.filename): self.pSV = polys
        self.theta = np.arctan(self.dz/np.sqrt(self.dx**2 + self.dy**2))
        self.output_dir = ''.join([DIR,'output/'])
        log.info(self.filename + " initialized")

     def find_flows():
        #Replace trunk with single node 

        #Remove edges with out flow 
        #Find the connected component with root in it, this is stem flow 

        #Remove all stem flow from orig graph, these are drip paths
            #Find node of minimal height in each component, this is the drip point 
    
    def network_simplex():
        #can be used to calculate flows on graphs with demands 
        #we could set a demand of generates X volume of flow 
