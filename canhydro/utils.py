from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd
import shapely.geometry as geometry
from scipy.spatial import Delaunay
from shapely.ops import polygonize

import matpotlib.pyplot as plt

from canhydro.global_vars import input_dir, log


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


def concave_hull(boundary_points, alpha):
    """alpha shape / concave hull
    Draws a minimal concave polygon with a concavity factor alpha"""
    if len(boundary_points) < 4:
        # When you have a triangle, there is no sense in computing an alpha
        # shape.
        return geometry.MultiPoint(list(boundary_points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):
        # adds a line between points i and j
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    coords = np.array([point.coords[0] for point in boundary_points])

    # Minimal set of triangles with points in set
    tri = Delaunay(coords)

    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]

        # Lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)

        # Semiperimeter of triangle
        s = (a + b + c) / 2.0

        # Area of triangle by Heron's formula
        area = math.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)

        # Here's the radius filter.
        # print circum_r
        if circum_r < 1.0 / alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)

    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return unary_union(triangles), edge_points
