from __future__ import annotations

import math
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely.geometry as geometry
from scipy.spatial import Delaunay,distance
from shapely.ops import polygonize, unary_union
from shapely.geometry import Polygon

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

    coords = np.array([point.coords[0] for point in boundary_points if len(point.coords)>0])

    # Minimal set of triangles with points in set
    tri = Delaunay(coords)

    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
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

def union(list):
    return unary_union(list)

def furthest_point(point, points):
    furthest_index = distance.cdist([point], points).argmax()
    return points[furthest_index]

def get_projection(vector:list, magnitude:list, radius:float()):
    noCirPoints = 360
    tCir = np.linspace(
        0, 2 * np.pi, noCirPoints
    )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
    a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
    b_ortho = np.sin(tCir)  # y coordinates of the points on a circle
    delt_a = magnitude[0]
    delt_b = magnitude[1]
    delt_c = magnitude[2]
    dim_a = vector[0]
    dim_b = vector[1]
    dim_c = vector[2]
    # unit vector at base of cylinder, pointing up cylinder axis
    vNorm = np.sqrt(delt_a**2 + delt_b**2 + delt_c**2)
    aV = np.hstack((delt_a, delt_b, delt_c)) / vNorm
    bV = -aV  # unit vector looking down from top circle (but not translated)
    # function to find orthgonal vectors
    oVz = lambda v, a, b: ((-v[0] * a - v[1] * b) / v[2])
    # initializing min max arrays+
    min_c = np.zeros_like(delt_c)
    max_c = np.zeros_like(delt_c)
    pSV = []
    projection = {
                    "polygon": Polygon(),
                    "base_vector": None,
                    "anti_vector": None,
                    "angle": None,
                    "area": None,
                }
    try:
        # for each cylinder
        if not np.isnan(dim_a[0]):
            if np.logical_and(delt_a == 0, delt_b == 0):
                pX = dim_a[0] + radius * a_ortho
                pY = dim_b[0] + radius * b_ortho
                cPS = Polygon(list(zip(pX, pY)))
                min_c = np.min(dim_c[:])
                max_c = np.max(dim_c[:])
            else:
                # find orthogonal vectors @ endpoints
                # Identifies corners of projected rectangle
                aVp1 = np.hstack((aV[1], -aV[0]))
                aVp2 = np.hstack((-aV[1], aV[0]))
                bVp1 = np.hstack((bV[1], -bV[0]))
                bVp2 = np.hstack((-bV[1], bV[0]))
                aVp1 = aVp1 / np.linalg.norm(aVp1)
                aVp2 = aVp2 / np.linalg.norm(aVp2)
                bVp1 = bVp1 / np.linalg.norm(bVp1)
                bVp2 = bVp2 / np.linalg.norm(bVp2)
                # from each endpoint, use radius to find vertices of the rectangle
                x1 = dim_a[0] + radius * aVp1[0]
                y1 = dim_b[0] + radius * aVp1[1]
                x2 = dim_a[0] + radius * aVp2[0]
                y2 = dim_b[0] + radius * aVp2[1]
                x3 = dim_a[1] + radius * bVp1[0]
                y3 = dim_b[1] + radius * bVp1[1]
                x4 = dim_a[1] + radius * bVp2[0]
                y4 = dim_b[1] + radius * bVp2[1]
                # calculate set of orthgonal vectors using lambda function
                # That is 360 orthogonal vectors ending at eqidistant points along
                # a circle of radius radius with the starting point of our cylinder
                # at is center
                ZOrtho = oVz(aV[:], a_ortho, b_ortho)
                # unit-ify the orthgonal vectors
                uovd = np.sqrt(a_ortho**2 + b_ortho**2 + ZOrtho**2)
                # Confounded - why does removing the first three [:,None]'s below lead to non-circular projections
                # for XZ?
                uov = (
                    np.hstack((a_ortho[:, None], b_ortho[:, None], ZOrtho[:, None]))
                    / uovd[:, None]
                )
                # donot re unit-fy, you only want the horizontal component, not the
                # renormalized horizontal component
                # using only the X and Y components, find circle coods in plane of
                # interest
                xaC = dim_a[0] + uov[:, 0] * radius
                yaC = dim_b[0] + uov[:, 1] * radius
                zaC = dim_c[0] + uov[:, 2] * radius
                xbC = dim_a[1] + uov[:, 0] * radius
                ybC = dim_b[1] + uov[:, 1] * radius
                zbC = dim_c[1] + uov[:, 2] * radius
                min_c = np.min(np.vstack((zaC, zbC)))
                max_c = np.max(np.vstack((zaC, zbC)))
                # assymble total package
                rX = np.vstack((x1, x2, x3, x4))
                rY = np.vstack((y1, y2, y3, y4))
                # test for circle parts in polygon
                try:
                    c1 = Polygon(
                        list(
                            zip(
                                [0 if np.isnan(x) else x for x in xaC],
                                [0 if np.isnan(y) else y for y in yaC],
                            )
                        )
                    )
                    bBox = Polygon(
                        list(
                            zip(
                                [0 if np.isnan(x) else x for x in rX],
                                [0 if np.isnan(y) else y for y in rY],
                            )
                        )
                    )
                    c2 = Polygon(
                        list(
                            zip(
                                [0 if np.isnan(x) else x for x in xbC],
                                [0 if np.isnan(y) else y for y in ybC],
                            )
                        )
                    )
                    partsPS = [c1, bBox, c2]
                except:
                    log.info("Error creating projection polygons")
                try:
                    cPS = unary_union(partsPS)
                except:
                    print(np.any(np.isnan(xaC)))
                    log.info("Error unioning projection polygons")
                    print(yaC)
                    print(rX)
                    print(rY)
                    print(xbC)
                    print(ybC)
                # get angle away from plane projected on to
                run = math.sqrt(delt_b**2 + delt_a**2)
                rise = delt_c
                if run == 0:
                    slope = 1  # straightDown e.g. is in flow
                else:
                    slope = rise / run
                ang = np.arctan(slope)
                projection = {
                    "polygon": cPS,
                    "base_vector": aV,
                    "anti_vector": bV,
                    "angle": ang,
                    "area": cPS.area,
                }
                return projection
        else:
            log.info(f"dim_a[0] is null, unable to project")
        return projection
    except UnboundLocalError:
        log.info(f"UnboundLocalError: vector : {vector} magnitude: {magnitude} radius: {radius}")


def countlines(start, lines=0, header=True, begin_start=None):
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