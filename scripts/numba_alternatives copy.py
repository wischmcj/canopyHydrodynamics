#utils 

import numpy as np
# from src.canhydro.geometry import numba_get_projection

from src.canhydro.CylinderCollection import CylinderCollection
from src.canhydro.Cylinder import Cylinder

from src.canhydro.utils import intermitent_log


from numba import njit, prange
from numba.typed import List



def stack(to_stack:list[np.array], col: bool = True):
    """
        A wrapper for njit stack that handles errors and allows for 
        less strict typing 
    """
    list_of_array = list(to_stack)
    try:
        njit_stack(list_of_array,col)
    except ValueError as err:
        left_shape = list_of_array[0].shape[0]
        right_shape = list_of_array[1].shape[0]
        stack_type ='column' if col else 'row'
        msg = f'{err}: Cannot {stack_type} stack arrays with shapes {left_shape} and {right_shape}'
        log.error(msg)
        raise ValueError(msg) from err

@njit()
def njit_stack(list_of_array:np.array[np.array()], col: bool):
    """
        numba doesn't play well with np stacks, so I had to do it myself
    """
    num_in = len(list_of_array)
    left_shape = list_of_array[0].shape[0]
    shape = (num_in, left_shape)
    stacked_array = np.empty(shape)
    for j in prange(len(list_of_array)): 
        stacked_array[j] = list_of_array[j]
    return stacked_array if not col else stacked_array.T

# geometry 

# @profile
@njit()

# *************** This is broken right now since casting a list of arrays as a numpy array is not no python friendly***************
# *****************https://stackoverflow.com/questions/53861099/creating-a-numpy-array-decorated-by-njit-from-numba
def projection_jit(vector: np.array, magnitude: np.array, radius: np.float32):
    dim_a = vector[0]
    dim_b = vector[1]
    dim_c = vector[2]
    delt_a = magnitude[0]
    delt_b = magnitude[1]
    delt_c = magnitude[2]

    noCirPoints = 360
    tCir = np.linspace(
        0, 2 * np.pi, noCirPoints
    )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
    a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
    b_ortho = np.sin(tCir)  # y coordinates of the points on a circle
    vNorm = np.sqrt(delt_a**2 + delt_b**2 + delt_c**2)
    # aV =  np.array([0.0,0.0,0.0]) if vNorm == 0 else np.array([delt_a, delt_b, delt_c])/ vNorm
    aV = np.array([delt_a, delt_b, delt_c]) / vNorm
    bV = -aV  # unit vector looking down from top circle (but not translated)
    # unit vector at base of cylinder, pointing up cylinder axis

    # function to find orthgonal vectors
    oVz = lambda v, a, b: ((-v[0] * a - v[1] * b) / v[2])
    # initializing min max arrays+
    pSV = []
    run = math.sqrt(delt_b**2 + delt_a**2)
    rise = delt_c
    if run == 0:
        slope = 1  # straightDown e.g. is in flow
    else:
        slope = rise / run
    ang = np.arctan(slope)
    print(f"projecting vector {aV}")
    # try:

    c1 = np.zeros((360, 2))
    # bBox = stack((
    #                     np.array([0.0,0.0,0.0,0.0]),
    #                     np.array([0.0,0.0,0.0,0.0])
    #                 ),
    #                 True
    #         )
    c2 = np.zeros((360, 2))
    bBox = np.zeros((4, 2))
    if not np.isnan(dim_a[0]):
        if np.logical_and(delt_a == 0.0, delt_b == 0.0):
            pX = dim_a[0] + radius * a_ortho
            pY = dim_b[0] + radius * b_ortho
            c2 = stack(np.array((pX, pY)), col=True)
            ang = 1.5708  # 90* straight up and down
        else:
            if aV[2] != 0.0:
                # calculate circular sections (only relevant if cyl is not parallel with the xy plane)
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

                # c1 =  np.column_stack((
                #            [0 if np.isnan(x) else x for x in xaC],
                #            [0 if np.isnan(y) else y for y in yaC],
                #         ))
                c1 = stack(
                    (
                        np.array([0 if np.isnan(x) else x for x in xaC]),
                        np.array([0 if np.isnan(y) else y for y in yaC]),
                    ),
                    col=True,
                )
                print(c1)
                xC = np.array([0 if np.isnan(x) else x for x in xbC])
                yC = np.array([0 if np.isnan(y) else y for y in ybC])
                # c2 =  np.column_stack((
                #             xC,
                #             yC,
                #         ))
                c2 = stack(np.array((xC, yC)), col=True)
            # calculating the rectangular portion of the projection
            # relevant for all cyls

            # find orthogonal vectors @ endpoints
            # Identifies corners of projected rectangle
            aVp1 = np.array([aV[1], -aV[0]])
            aVp2 = np.array([-aV[1], aV[0]])
            bVp1 = np.array([bV[1], -bV[0]])
            bVp2 = np.array([-bV[1], bV[0]])
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

            # bBox =  np.array(
            #         [   np.array([x1,y1]),
            #             np.array([x2,y2]),
            #             np.array([x3,y3]),
            #             np.array([x4,y4])
            #         ]
            #         )
            # bBox =  np.array(
            #         (  (x1,y1),
            #            (x2,y2),
            #            (x3,y3),
            #            (x4,y4)
            #         )
            #         )
            bBox = stack(np.array((np.array([x1, x2, x3, x4]), np.array([y1, y2, y3, y4]))), True)

            # breakpoint()
        # print(".")
        # print(typeof(bBox))
        # print(".")
        # print(typeof(c1))
        # print(".")
        # print(typeof(c2))
        # print(".")
        # print(typeof(aV))
        return c1, bBox, c2, ang, aV, bV
    # else:
    #     return None
    # except UnboundLocalError:
    #     log.info(
    #         f"UnboundLocalError: vector : {vector} magnitude: {magnitude} radius: {radius}"
    #     )

@njit()
def numba_get_projection(vector: list, magnitude: list, radius: np.float32):
    """
    Takes in the vector (starting point), magnitude and radius that fully define a cylinder.
    Finds the projection of the cylinder on a plane
    """
    projection = {
        "polygon": Polygon(),
        "base_vector": (0, 0, 0),
        "anti_vector": (0, 0, 0),
        "angle": 0,
        "area": 0,
    }
    c1, bBox, c2, ang, aV, bV = projection_jit(
        np.array(vector), np.array(magnitude), np.float32(radius)
    )
    partsPS = [c1, bBox, c2]
    if np.max([poly_part.size for poly_part in partsPS]) > 0:
        to_union = [arr for arr in partsPS if arr.size > 0]
        cPS = unary_union([Polygon(part) for part in to_union])
        projection = {
            "polygon": cPS,
            "base_vector": aV,
            "anti_vector": bV,
            "angle": ang,
            "area": cPS.area,
        }

    # if c1.size > 0:
    #     try:
    #         c1 = Polygon(c1)
    #         bBox = Polygon(bBox)
    #         c2 = Polygon(c2)
    #         partsPS = [c1, bBox, c2]
    #     except:
    #         log.info("Error creating projection polygons")
    #     try:
    #         cPS = unary_union(partsPS)
    #     except:
    #         log.info("Error unioning projection polygons")
    #     # get angle away from plane projected on to
    # else:
    # cPS = Polygon(c2)

    return projection

# Cylinder Collection 
class NumbaCylinderCollection(CylinderCollection):

    def __init__(self, **kwargs):
        super().__iniit__(**kwargs)    

    def numba_project_cylinders(self, plane: str = "XY", force_rerun: bool = False):
        """
        This a wrapper to the workhorse function:
        numba_get_projection. Numba is a tool that compiles 
        python code into machine code, which can be much faster.
        
        This function is still under development, however, and
        appears to have minimal/no advantage over the 
        non-numba version
        """
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        elif not force_rerun and self.projections[plane]:
            log.info(
                "cached pro     jections exist, pass 'force_rerun=True to calculate new projections "
            )
        else:
            polys = []
            log.info(f"Projection into {plane} axis begun for file {self.file_name}")
            for idx, cyl in enumerate(self.cylinders):
                poly = cyl.numba_get_projection(plane)
                polys.append(poly)
                # print a progress update once every 10 thousand or so cylinders
                intermitent_log(idx, self.no_cylinders, "Cylinder projection: ")
            # Used by other functions to know what projections have been run
            self.projections[plane] = True
            self.pSV = polys

#Cylinder 
class NumbaCylinder(Cylinder):

    def __init__(self, **kwargs):
        super().__iniit__(**kwargs)    

    def numba_get_projection(self, plane="XY"):
        noCirPoints = 360
        tCir = np.linspace(
            0, 2 * np.pi, noCirPoints
        )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
        a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
        b_ortho = np.sin(tCir)  # y coordinates of the points on a circle
        if plane == "XY":
            magnitude = [self.dx, self.dy, self.dz]
            vector = [np.transpose(self.x), np.transpose(self.y), np.transpose(self.z)]
        elif plane == "XZ":
            magnitude = [self.dx, self.dz, self.dy]
            vector = [np.transpose(self.x), np.transpose(self.z), np.transpose(self.y)]
        else:
            magnitude = [self.dy, self.dz, self.dx]
            vector = [np.transpose(self.y), np.transpose(self.z), np.transpose(self.x)]
        projection = numba_get_projection(vector, magnitude, self.radius)
        self.projected_data[plane] = projection
        if plane == "XY":
            self.xy_area = self.projected_data["XY"]["area"]
        return projection["polygon"]


#efficiency_test

# #@profile
# def test_small_tree_proj(small_tree, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     small_tree.project_cylinders("XY")
#     assert 1 == 1


# #@profile
# def test_numba_small_tree_proj(small_tree, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     small_tree.numba_project_cylinders("XY")
#     assert 1 == 1


# #@profile
# def test_happy_proj(happy_path_projection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     happy_path_projection.project_cylinders("XZ")
#     assert 1 == 1


# #@profile
# def test_numba_happy_proj(happy_path_projection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     happy_path_projection.numba_project_cylinders("XZ")
#     assert 1 == 1


# #@profile
# def test_large_proj(large_collection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     large_collection.project_cylinders("XZ")
#     assert 1 == 1


# #@profile
# def test_numba_large_proj(large_collection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     large_collection.numba_project_cylinders("XZ")
#     assert 1 == 1
