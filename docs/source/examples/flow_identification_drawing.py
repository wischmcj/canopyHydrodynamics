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

# %% [markdown]
# # Use Case Examples
#     This notebool is intended to provide snippets of code for various use cases that users might have for the canoPyHydro package
# %%
# Reading in a QSM and projecting its cylinders from 3D to 2D
from canoPyHydro.CylinderCollection import CylinderCollection

# Initializing a CylinderCollection object
myCollection = CylinderCollection()

# Converting a specified file to a CylinderCollection object
myCollection.from_csv("example_tree.csv")

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.project_cylinders("XZ")
myCollection.draw("XZ")

# Requesting an plot of the tree projected onto the YZ plane ('side' view)
myCollection.project_cylinders("YZ")
myCollection.draw("YZ")

# Requesting an plot of the tree projected onto the XY plane (Birds eye view)
myCollection.project_cylinders("XY")
myCollection.draw("XY")

# %%
# Finding and Highlighting the stem flow component

import os

os.environ["CANOPYHYDRO_CONFIG"] = "./canopyhydro_config.toml"
from src.canopyhydro.CylinderCollection import CylinderCollection

# Initializing a CylinderCollection object
myCollection = CylinderCollection()

# Converting a specified file to a CylinderCollection object
myCollection.from_csv("charlie_brown.csv")

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.project_cylinders("XY")

# creating the digraph model
myCollection.initialize_digraph_from()
# Identifying the flows to which each cyl belongs
myCollection.find_flow_components()
# Calculating the propreties of each flow
myCollection.calculate_flows()


myCollection.draw(
    "XZ", highlight_lambda=lambda: is_stem, save=True, file_name_ext="docs_ex"
)

# %%
# Drawing alpha shapes/convex hulls around the tree

from src.canopyhydro.CylinderCollection import CylinderCollection

# Creating a CylinderCollection object
myCollection = CylinderCollection()
myCollection.from_csv("5_SmallTree.csv")
myCollection.project_cylinders("XY")
myCollection.initialize_digraph_from()
myCollection.find_flow_components()
myCollection.calculate_flows()

# drawing the tree for reference
myCollection.draw("XY", save=True, file_name_ext="read_me_alpha")

# Drawing the whole canopy boundary
myCollection.watershed_boundary(plane="XY", draw=True)

# Drawing the whole canopy boundary with a looser fit
myCollection.watershed_boundary(
    plane="XY",
    curvature_alpha=0.3,  # determines how 'tight' the fit is
)

# Saving the whole canopy boundary plot to a file
myCollection.watershed_boundary(
    plane="XY", save=True, file_name_ext="read_me_alpha", draw=True
)
# Plotting the most recent canopy boundary as an overlay on the tree's projected cylinders
myCollection.draw(plane="XY", include_alpha_shape=True)

# %%
# Plotting projections; highlighting and filtering
# Plot the entire tree with stem flow highlighted
myCollection.draw("XZ", highlight_lambda=lambda: is_stem)

# Plot the interesting portion of the tree with stem flow highlighted
myCollection.draw(
    "XZ", highlight_lambda=lambda: is_stem, filter_lambda=lambda: cyl_id > 100
)

# Adding drip points to the above mentioned plot
myCollection.draw(
    "XZ",
    highlight_lambda=lambda: is_stem,
    filter_lambda=lambda: cyl_id > 100,
    include_drips=True,
)

# %%
# Extracting and printing flow data

from tomark import Tomark


def round_if(to_round, n):
    if type(to_round) == tuple:
        return tuple(round(x, 1) for x in to_round)
    else:
        return round(to_round, n)


data = [f.__dict__ for f in myCollection.flows]
data = [d for d in data if d["num_cylinders"] > 10]
data = [{k: round_if(v, 3) for k, v in d.items()} for d in data]
markdown = Tomark.table(data)
print(markdown)

# %%
# Plotting a trees cylinders, watershed boundary and drip points all together

myCollection.watershed_boundary(
    plane="XY",
    curvature_alpha=0.15,
    filter_lambda=lambda: is_stem and cyl_id > 100,
    draw=True,
)
# Adding drip points to the above mentioned plot
myCollection.draw(
    "XY",
    include_alpha_shape=True,
    highlight_lambda=lambda: is_stem,
    filter_lambda=lambda: cyl_id > 100,
)

# %%
# Drawing subsets of a tree's cylinders based on different QSM data

import os

from matplotlib import pyplot as plt

os.environ["CANOPYHYDRO_CONFIG"] = "./canopyhydro_config.toml"
from src.canopyhydro.CylinderCollection import CylinderCollection

# Initializing a CylinderCollection object
myCollection = CylinderCollection()

# Converting a specified file to a Cyl0nderCollection object
myCollection.from_csv("example_tree.csv")

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.project_cylinders("XY")

ax[1] = myCollection.draw(
    "XZ",
    filter_lambda=lambda: branch_order <= 2,
    highlight_lambda=lambda: branch_order == 2,
    save=True,
    file_name_ext="_bo_le_2",
)
ax[2] = myCollection.draw(
    "XZ",
    filter_lambda=lambda: branch_order <= 3,
    highlight_lambda=lambda: branch_order == 3,
    save=True,
    file_name_ext="_bo_le_3",
)
ax[3] = myCollection.draw(
    "XZ",
    filter_lambda=lambda: branch_order <= 4,
    highlight_lambda=lambda: branch_order == 4,
    save=True,
    file_name_ext="_bo_le_4",
)
plt.show()
