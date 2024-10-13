# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: venv
#     language: python
#     name: python3
# ---

# %%
from __future__ import annotations

import functools

import mpl_toolkits.mplot3d.art3d as art3d
# Getting rotation matrix
import numpy as np
from matplotlib.patches import Circle

# %% [markdown]
# The below is intended to help users visualize how canopyHydro transforms 3D cylinders (from QSM models) to 2D shapes.
# %%
# The functions in this block are used to define a 3D cylinder in space





def rotate_z_to_vector(b: np.array):
    if np.linalg.norm(b) == 0:
        return np.eye(3)
    if np.linalg.norm(b) != 1:
        raise ValueError("b must be a unit vector")
    # Algorithm from https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation
    # b must be unit vector
    # a is the z unit vector
    a = [0, 0, 1]
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)
    # The skew-symmetric cross product matrix of v
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) * -1
    # Rotation matrix as per Rodregues formula
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s**2))
    return R


@functools.lru_cache
def get_unoriented_cylinder(r, h, a=0, noCirPoints=200, nv=200):
    """
    Returns the parameterization of a cylinder given the radius (r), height (h), and origin (a)
    """
    theta = np.linspace(0, 2 * np.pi, noCirPoints)
    v = np.linspace(a, a + h, nv)
    theta, v = np.meshgrid(theta, v)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = v
    return x, y, z


def get_cylinder(r, h, a=0, noCirPoints=200, nv=200):
    """
    Returns the parameterization of a cylinder given the radius (r), height (h), and origin (a)
    """
    theta = np.linspace(0, 2 * np.pi, noCirPoints)
    v = np.linspace(a, a + h, nv)
    theta, v = np.meshgrid(theta, v)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = v
    #     rotation_matrix = np.array([[cos(a), -sin(a)], [sin(a), cos(a)]])
    # x   , y = zip(*[(x,y) @ rotation_matrix for x,y in zip(x,y)])

    return x, y, z


# %%
# Here we use the above functions to draw a 3d Cylinder aligned with z axis
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

Xc, Yc, Zc = get_cylinder(2, 4)
# end_cap = [Xc[0], Yc[0], Zc[0]]

ax.plot_surface(Xc, Yc, Zc, alpha=0.5)
plt.show()

import functools

import mpl_toolkits.mplot3d.art3d as art3d
# %%
# Below we draw a 3D cylinder aligned with an arbitrary vector
#  as well as the 2D projections of the cylinder on the xy, xz, and yz planes
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

v = np.array([3, 1, 7])

vu = v / np.linalg.norm(v)

Xc, Yc, Zc = get_cylinder(2, 4)

R = rotate_z_to_vector(vu)
t = np.transpose(np.array([Xc, Yc, Zc]))
x, y, z = np.transpose(t @ R, (2, 0, 1))

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
zmin = -4

ax.plot_surface(x, y, z, alpha=0.7)

# Drawing the vector v
ax.quiver(0, 0, 0, 1.5, 0.5, 3.5, color="black")

# ax.contourf(x, y, z, zdir='z', offset=-1, cmap='Greys', alpha = .5) #'coolwarm')
ax.contourf(x, y, z, zdir="z", offset=zmin, colors="C0")

# Drawing a circle to fill in gap left by contourf
# represents the end cap of the cylinder
# hard coded for this example for simplicity
p = Circle((0.5, 0.15), 1.77)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=zmin, zdir="z")  # , color = 'Grey')


plt.contourf(x, y, z, zdir="x", offset=zmin, colors="C0", labels="YZ")
# yz_contour.collections[0].set_label('YZ')
ax.contourf(x, y, z, zdir="y", offset=6, colors="C0")


ax.view_init(elev=30, azim=-45, roll=0)
ax.set(xlim=(-4, 6), ylim=(-5, 6), zlim=(zmin, 6), xlabel="X", ylabel="Y", zlabel="Z")

print(dir(ax))
plt.show()

# %%
