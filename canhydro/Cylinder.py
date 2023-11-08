"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import calendar
import os
from multiprocessing import Pool
from pathlib import Path
from random import random
from time import sleep
from math import pi, sqrt, isnan

from canhydro.global_vars import (
    qsm_cols,
    input_dir,
    output_dir,
    log
)

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import toml
from matplotlib.pyplot import cm

from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
# from descartes import PolygonPatch
# from mpl_toolkits import mplot3d


# import time
# import copy
# import math
# import openpyxl
# import geopandas as geo
NAME = "Cylinder"


class Cylinder:
    # A cylinder could arguably be a data class, it also could extend Shapely.polygon
    #
    # Its purpose is to act as a 1-1 copy of the rows in the QSM data
    # (rows that represent cylinders in 3-D space)
    #
    #

    # initialize our object level variables for cylider objects
    def __init__(self) -> None:
        # Base Attributes from file read in
        # self._cylinderCollection = CylinderCollection
        self.x = np.nan  # len 2 array
        self.y = np.nan  # len 2 array
        self.z = np.nan  # len 2 array
        self.dx = np.nan  # len 2 array
        self.dy = np.nan  # len 2 array
        self.dz = np.nan  # len 2 array
        self.radius = np.nan
        self.length = np.nan
        self.branch_order = np.nan
        self.branch_id = np.nan
        self.volume = np.nan
        self.parent_id = np.nan
        self.rev_branch_order = np.nan
        self.segment_id = np.nan

        self.surface_area = np.nan
        self.sa_to_vol = np.nan
        self.slope = np.nan

        # Calculated based off of projection
        self.projected_data = {
            "XY": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            },
            "YZ": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            },
            "XZ": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            }
        }

        # graph attributes, all in 3-D
        self.flow_id = np.nan
        self.flow_type = np.nan
        self.begins_at_drip_point = np.nan  # bool
        self.begins_at_divide_point = np.nan  # bool

        # others
        self.stem_path_id = np.nan

    def calc_surface_area(self):
        radius = self.radius
        length = self.length
        sa = 2*np.pi*radius*(radius + length) - 2*np.pi*radius*radius
        return sa

    def to_dict(self, cols:str = ''):
        if cols == '':
            """creates a cylinder corrosponding to that defined by a given row of the qsm (attrs)"""
            attr_dict = {
                    "x":self.x  
                    ,"y":self.y  
                    ,"z":self.z  
                    ,"dx":self.dx 
                    ,"dy":self.dy 
                    ,"dz":self.dz 
                    ,"radius"         :self.radius             
                    ,"length"          :self.length             
                    ,"surface_area"     :self.surface_area       
                    ,"volume"           :self.volume             
                    ,"sa_to_vol"        :self.sa_to_vol          
                    ,"branch_order"     :self.branch_order       
                    ,"branch_id"        :self.branch_id          
                    ,"parent_id"       :self.parent_id          
                    ,"rev_branch_order" :self.rev_branch_order   
                    ,"segment_id"       :self.segment_id         
                    ,"angle"            :self.angle  
                }            
        return attr_dict

    def create_from_list(self, attrs:list, columns = qsm_cols):
        """creates a cylinder corrosponding to that defined by a given row of the qsm (attrs)"""

        extract = lambda attr : attrs[columns[attr]] #pulls a column from the qsm row (attrs) corrosponding to the input attribute

        self.x = [extract('x')[0], extract('x')[1]]
        self.y = [extract('y')[0], extract('y')[1]]
        self.z = [extract('z')[0], extract('z')[1]]
        self.dx = self.x[1]-self.x[0]
        self.dy = self.y[1]-self.y[0]
        self.dz = self.z[1]=self.z[0]
        self.radius = extract('radius')
        self.length =  extract('length')
        self.surface_area = self.calc_surface_area()
        self.volume = extract('volume')
        self.sa_to_vol = self.surface_area/self.volume 
        self.branch_order =  extract('branch_order')
        self.branch_id =  extract('branch_id')
        self.parent_id = extract('parent_id')
        self.rev_branch_order = extract('reverse_branch_order')
        self.segment_id = extract('segment_id')
        self.angle = np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
        log.info(str(self.to_dict()))

    def get_projection(self, plane= 'XY'):
        noCirPoints = 360
        tCir = np.linspace(0,2*np.pi,noCirPoints) #360 evenly spaced points between 0 - 2pi (radian degrees)
        a_ortho = np.cos(tCir) #x coordinates of the points on a circle
        b_ortho = np.sin(tCir) #y coordinates of the points on a circle      
        
        if plane == 'XY':
            delt_a = self.dx 
            delt_b = self.dy
            delt_c = self.dz
            dim_a = np.transpose(self.x)
            dim_b = np.transpose(self.y) 
            dim_c = np.transpose(self.z) 
        elif plane == 'XZ':
            delt_a = self.dx 
            delt_b = self.dz
            delt_c = self.dy
            dim_a = np.transpose(self.x)   
            dim_b = np.transpose(self.z) 
            dim_c = np.transpose(self.y)
        else:
            delt_a = self.dz 
            delt_b = self.dy
            delt_c = self.dx
            dim_a = np.transpose(self.z)   
            dim_b = np.transpose(self.y) 
            dim_c = np.transpose(self.x)

        #unit vector at base of cylinder, pointing up cylinder axis
        vNorm = (np.sqrt(delt_a**2+delt_b**2+delt_c**2))
        aV = np.hstack((delt_a[:,None],delt_b[:,None],delt_c[:,None]))/vNorm[:,None]
        bV = -aV; #unit vector looking down from top circle (but not translated)
        
        #function to find orthgonal vectors
        oVz = lambda v,a,b : ((-v[0]*a - v[1]*b)/v[2])

        #initializing min max arrays+
        min_c = np.zeros_like(delt_c)
        max_c = np.zeros_like(delt_c)

        pSV = []

        #for each cylinder
        if not np.isnan(dim_a[0]):
            if np.logical_and(delt_a == 0, delt_b == 0):
        
                pX = dim_a[0] + self.radius*a_ortho
                pY = dim_b[0] + self.radius*b_ortho
                cPS = Polygon(list(zip(pX,pY)))
                min_c = np.min(dim_c[:])
                max_c = np.max(dim_c[:])
            else:
                #find orthogonal vectors @ endpoints
                #Identifies corners of projected rectangle 
                aVp1 = np.hstack((aV[1],-aV[0]))
                aVp2 = np.hstack((-aV[1],aV[0]))
                bVp1 = np.hstack((bV[1],-bV[0]))
                bVp2 = np.hstack((-bV[1],bV[0]))
            
                aVp1 = aVp1/np.linalg.norm(aVp1)
                aVp2 = aVp2/np.linalg.norm(aVp2)
                bVp1 = bVp1/np.linalg.norm(bVp1)
                bVp2 = bVp2/np.linalg.norm(bVp2)
            
                #from each endpoint, use radius to find vertices of the rectangle
                x1 = dim_a[0] + self.radius*aVp1[0]
                y1 = dim_b[0] + self.radius*aVp1[1]
                x2 = dim_a[0] + self.radius*aVp2[0]
                y2 = dim_b[0] + self.radius*aVp2[1]  
                x3 = dim_a[1] + self.radius*bVp1[0]
                y3 = dim_b[1] + self.radius*bVp1[1]
                x4 = dim_a[1] + self.radius*bVp2[0]
                y4 = dim_b[1] + self.radius*bVp2[1]
        
                #calculate set of orthgonal vectors using lambda function
                ZOrtho = oVz(aV[:],a_ortho,b_ortho)
            
                #unit-ify the orthgonal vectors
                uovd = np.sqrt(a_ortho**2 + b_ortho**2 + ZOrtho**2)
                uov = np.hstack((a_ortho[:,None],b_ortho[:,None],ZOrtho[:,None]))/uovd[:,None]
        
                #donot re unit-fy, you only want the horizontal component, not the
                #renormalized horizontal component
            
                #using only the X and Y components, find circle coods in plane of
                #interest
                xaC = dim_a[0] + uov[:,0]*self.radius
                yaC = dim_b[0] + uov[:,1]*self.radius
                zaC = dim_c[0] + uov[:,2]*self.radius
        
                xbC = dim_a[1] + uov[:,0]*self.radius
                ybC = dim_b[1] + uov[:,1]*self.radius
                zbC = dim_c[1] + uov[:,2]*self.radius
        
                min_c = np.min(np.vstack((zaC, zbC)))
                max_c = np.max(np.vstack((zaC, zbC)))
            
                #assymble total package
                rX = np.vstack((x1, x2, x3, x4))
                rY = np.vstack((y1, y2, y3, y4))
        
                #test for circle parts in polygon
                try:
                    c1 = Polygon(list(zip([0 if isnan(x) else x for x in xaC],[0 if isnan(y) else y for y in yaC])))
                    bBox = Polygon(list(zip([0 if isnan(x) else x for x in rX],[0 if isnan(y) else y for y in rY])))
                    c2 = Polygon(list(zip([0 if isnan(x) else x for x in xbC],[0 if isnan(y) else y for y in ybC])))
                            
                    partsPS = [c1,bBox,c2]
                except:
                    log.info(f"Error creating projection polygons")
            
                try: 
                    cPS = unary_union(partsPS)
                    breakpoint()
                except:
                    print(np.any(np.isnan(xaC)))
                    log.info(f"Error unioning projection polygons")
                    print(yaC)
                    print(rX)
                    print(rY)
                    print(xbC)
                    print(ybC)
        breakpoint()
        # Calculated based off of projection
        self.projected_data[plane]['unit_vect'] = aV
        self.projected_data[plane]['rev_unit_vect'] = bV
        self.projected_data[plane]['rise_angle_rad'] = np.nan #havent thought through yet 
        self.projected_data[plane]['poly'] = cPS

        return pSV


    def get_flow_data():
        """Returns the flow ID and flow characteristics of the flow the cyl is contained in"""
        print("Get flow data not written")
