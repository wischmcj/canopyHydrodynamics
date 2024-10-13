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

# %%
#   RUNME.setup
#  This helps set up your kernal's environment in order to avoid errors
import os

import numpy as np
from matplotlib import pyplot as plt

from canopyhydro.Cylinder import Cylinder
from src.canopyhydro.CylinderCollection import CylinderCollection
from src.canopyhydro.Forester import Forester

# Determines where configuration file is located
# file contains directory info and model input settings
config_file = os.environ[
    "CANOPYHYDRO_CONFIG"
] = f"{os.getcwd()}/canopyhydro_config.toml"
log_config = os.environ["CANOPYHYDRO_LOG_CONFIG"] = f"{os.getcwd()}/logging_config.yml"


# %%
# Creating a Cylinder object via inputs
myCyl = Cylinder(
    cyl_id=1.0,
    x=[0, 3],
    y=[0, 2],
    z=[0, 6],
    radius=2.0,
    length=0.064433,
    branch_order=0.0,
    branch_id=0.0,
    volume=0.010021,
    parent_id=0.0,
    reverse_branch_order=32.0,
    segment_id=0.0,
)
print(myCyl)

# %% [markdown]
# ### Cylinder

# %% [markdown]
# The Cylinder class is used to represent the 3-D cylinders that make up a QSM

# %%
# A trivial example of a Cylinder object


myCyl = Cylinder(
    cyl_id=1.0,
    x=[3, 6],
    y=[2, 4],
    z=[6, 12],
    radius=2.0,
    length=0.064433,
    branch_order=0.0,
    branch_id=0.0,
    volume=0.010021,
    parent_id=0.0,
    reverse_branch_order=32.0,
    segment_id=0.0,
)

fig = myCyl.draw_3D(show=True, draw_vectors=True, draw_projections=True)

# %% [markdown]
# The most important function of Cylinder objects is their ability to return data regarding the projections onto planes. Cylinder objects utilize our custom 'geometry' module to calculate their projections onto the XY, XZ and YZ planes.

# %%
# Projecting a cylinder onto the XY plane
import numpy as np

from canopyhydro.Cylinder import Cylinder
from canopyhydro.CylinderCollection import CylinderCollection

cyl = Cylinder(
    1, np.array([0, 1]), np.array([0, 1]), np.array([0, 1]), 1, 1, 0, 0, 1, 0, 0, 0
)
cyl.get_projection("XY")
print(cyl.projected_data["XY"]["polygon"])
print(cyl.projected_data["XY"]["base_vector"])
print(cyl.projected_data["XY"]["anti_vector"])
print(cyl.projected_data["XY"]["angle"])
print(cyl.projected_data["XY"]["area"])

# %% [markdown]
# The main use of this functionality is shown below. Namely, projecting cylinders allows us to more readily provide visualizations of the tree canopy they represent

# %%
# Here we show the 3D view and the 3 possible 2D projections of a Cylinder object
myCyl = Cylinder(
    cyl_id=1.0,
    x=[0, 3],
    y=[0, 2],
    z=[0, 6],
    radius=2.0,
    length=0.064433,
    branch_order=0.0,
    branch_id=0.0,
    volume=0.010021,
    parent_id=0.0,
    reverse_branch_order=32.0,
    segment_id=0.0,
)


fig = myCyl.draw_3D(show=False, draw_projections=True)

myCyl.get_projection("XY")
print("'myCyl' as seen from above")
print(
    f"The 'XY' projection of myCyl as an area of {round(myCyl.projected_data['XY']['area'],2)} cm^3"
)
print(
    f"       and the cylinder makes an angle of {round(myCyl.projected_data['XY']['angle'],2)} radians with the XY plane"
)

myCyl.get_projection("XZ")
print(
    f"The 'XZ' projection of myCyl as an area of {round(myCyl.projected_data['XZ']['area'],2)} cm^3"
)
print(
    f"       and the cylinder makes an angle of {round(myCyl.projected_data['XZ']['angle'],2)} radians with the XZ plane"
)

myCyl.get_projection("YZ")
print(
    f"The 'YZ' projection of myCyl as an area of {round(myCyl.projected_data['YZ']['area'],2)} cm^3"
)
print(
    f"       and the cylinder makes an angle of {round(myCyl.projected_data['YZ']['angle'],2)} radians with the YZ plane"
)

# %% [markdown]
# In the above 3D representations, the entire surface of cylinder must be calculated, s this may be computationally intensive. \
# As such, the 'get_projection' function calculates statistics regarding 2D projections directly instead.

# %%
# The get_projection function allows for the retrieval of
# projection data without the need to define the entire surface of the cylinder


myCyl = Cylinder(
    cyl_id=1.0,
    x=[3, 6],
    y=[2, 4],
    z=[6, 12],
    radius=2.0,
    length=0.064433,
    branch_order=0.0,
    branch_id=0.0,
    volume=0.010021,
    parent_id=0.0,
    reverse_branch_order=32.0,
    segment_id=0.0,
)

print("'myCyl' as seen from above")
myCyl.get_projection("XY")
myCyl.draw(plane="XY")
plt.show()

print("'myCyl' as seen from the 'fromt' of the tree")
myCyl.get_projection("XZ")
myCyl.draw(plane="XZ")
plt.show()

print("'myCyl' as seen from one 'side' of the tree")
myCyl.get_projection("YZ")
myCyl.draw(plane="YZ")
plt.show()

# %% [markdown]
# For more information regarding these 2D projections, see [Projecting Cylinders](projecting_cylinders.ipynb)

# %% [markdown]
# ### Cylinder Collection

# %% [markdown]
# Cylinder Collections are just as they sound and, at the most basic level, a Cylinder Collection is defined as a list of 1 or more Cylinder objects. \
# Cylinder Collections almost always represent QSM's (or parts of a QSM), and are meant to help users explore their QSMs. \
# Below, we demonstrate how one might initialize a cylinder collection using cylinder data (e.g. QSM data) stored in a CSV file.

# %%
# Example showing the most basic possible cylinder Collection4
myCollection = CylinderCollection()
# The below file is one of our several testing files, featuring only
# the trunk of a tree and one of its branches
myCollection.from_csv("5_SmallTree.csv")
myCollection.statistics("XY")
myCollection.statistics("XY")

# %% [markdown]
#
# The below demonstrates how these projection and drawing functionality of the cylinder class can be extended to allow for drawing, highlighting and filtering entire tree canopies
#

# %%
myCollection = CylinderCollection()
# The below file is one of our several testing files, featuring only
# the trunk of a tree and one of its branches
myCollection.from_csv("10_MediumCollection.csv")

# Collections can be drawn without any frills
myCollection.draw("XY", show=True)

# One can also choose to draw only a portion of the collection,
# utilizing 'lambda' functions to define the desired portion of the tree
myCollection.draw("XZ", show=False, filter_lambda=lambda: branch_id < 20)


# Similarly, features of interest can be highlighted in a different color
myCollection.draw("XZ", show=False, filter_lambda=lambda: cyl_id > 100)

myCollection.draw(
    "XZ",
    filter_lambda=lambda: cyl_id > 50,
    highlight_lambda=lambda: branch_order > 0,
    save=True,
    file_name_ext="highlighted_branch_tutorial.svg",
)  # noqa


myCollection.draw(
    "XZ",
    filter_lambda=lambda: cyl_id > 50,
    highlight_lambda=lambda: cyl_id > 100,
    save=True,
    file_name_ext="highlighted_branch_tutorial.svg",
)


# %% [markdown]
# We recommend sticing to 2D drawings where possible, as this is a lot less computationally intensive than the alternative. However the same '.draw' function can even generate 3D plots of the tree!

# %%
myCollection = CylinderCollection()
# The below file is one of our several testing files, featuring only
# the trunk of a tree and one of its branches
myCollection.from_csv("charlie_brown.csv")

# by filtering for cyl_id>100, we are only plotting the
# cylinders that are part of the branch
# myCollection.draw("XZ", show=True, save=True)
# myCollection.draw("XY", show=True, save=True)
# myCollection.draw("YZ", show=True, save=True)
myCollection.draw("3D", show=False, save=True, file_name_ext="2d_3d_comparison")

print("XZ Projection of a collection of cylinders")

# %% [markdown]
# To understand how these features support canoPyHydro's titular 'flow identification' functionality, check out [Flow Identification and Drawing](flow_identification_drawing.ipynb)

# %% [markdown]
# #

# %% [markdown]
# ### Forester

# %% [markdown]
# Forester objects allow users to conveniently create and manage Cylinder Collections. In particular, Foresters are useful for reading in and processing QSM files.
#

# %% [markdown]
# When a Forester object is created, available file names are read from the default directory, './data/input/'. \
# This list of available files can be accessed through the Forester.file_names attribute, as shown below

# %% tags=["myTag"]
# Creating a new Forester object
myForester = Forester()
print(
    f"Files available in {myForester.directory}: {list(map(str,myForester.file_names))}"
)

# %% [markdown]
# Optionally, a custom path may be passed to the Forester object, In which case, the Forester will look for files in the passed directory instead

# %%
# Passing a custom directory to the Forester object will change the directory attribute
directory = "/data/test/"
myForester = Forester("data/test/")
print(
    f"Files available in {myForester.directory}: {list(map(str,myForester.file_names))}"
)

# %% [markdown]
# The 'qsm_to_collection' function can be used create CylinderCollections from a specified file.

# %%
# Importing a QSM file as a CylinderCollection
myForester = Forester("data/test/")
myForester.qsm_to_collection("5_SmallTree.csv")

cylCollections = myForester.cylinder_collections
firstCollection = cylCollections[0]

print(
    f"Forester has {len(cylCollections)} CylinderCollection, imported from {cylCollections[0].file_name}"
)

# %% [markdown]
# If 'All' is provided as the file name, all of the files in the given directory will be read in as CylinderCollections. \
# (Note that this may require a significant amount of memory.)

# %%
# Reading in all files in the directory as collections
myForester.qsm_to_collection("All")
cylCollections = myForester.cylinder_collections
firstCollection = cylCollections[0]
print(
    f"""Forester created {len(cylCollections)} CylinderCollections, imported from the following files
      {list(map(lambda x: x.file_name,cylCollections))}"""
)

# %%
# Code stored here that may or may not be useful as scraps
# Alternatively, the 'Forester' class can be used
myForester = Forester()
print(
    f"Files available in {myForester.directory}: {list(map(str,myForester.file_names))}"
)
# Foresters can ...
## ... use a custom directory
myForester = Forester("data/test/")
print(
    f"Files available in {myForester.directory}: {list(map(str,myForester.file_names))}"
)
## ... Read in single QSMs
myForester.qsm_to_collection("5_SmallTree.csv")
print(f"Forester has {len(myForester.cylinder_collections)} CylinderCollections")
# .. or read in all files in a directory
myForester.qsm_to_collection("All")
print(
    f"""Forester created {len(myForester.cylinder_collections)} CylinderCollections,"""
)

# %% [markdown]
# Putting all that we have learned together, you can see that all of the statistics available through canoPyHydro can be generated with the below 10 lines of code

# %%
forest = Forester(test_input_dir)
forest.get_file_names()
forest.qsm_to_collection(file_name="3_HappyPathWTrunk.csv")
collection = forest.cylinder_collections[0]
collection.project_cylinders("XY")
collection.initialize_digraph_from(in_flow_grade_lim=-0.16)
collection.find_flow_components()
print("finished_find_flow_components")
collection.calculate_flows()
pickle_collection(collection)
