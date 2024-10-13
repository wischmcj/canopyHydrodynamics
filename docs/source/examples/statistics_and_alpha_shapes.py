# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.4
# ---

from __future__ import annotations

# %% [markdown]
#
# The most important of these attributes is CylinderCollection.digraph, which is a mathematical graph corresponding the CylinderCollection. \
# This graph representation is used in tandem with a traversal algorithm to predict which rain partition each cylinder in the collection belongs
# %%

# Code example printing out a cylinder collection, colored by stem v. drip flow


# %% [markdown]
#
# Alpha shapes are another key attribute used in statistics calculations. \
# Alpha shapes represent the estimated area covered by the represented tree's canopy when projected in the XY, XZ or YZ direction. \
# **these shapes are particularly important in the calculation of Woody Index (see statistics_calculations for more info)
#

# %%
# A demonstration showing the calculation and plotting of alpha shapes


# %% [markdown]
#
# The remaining attributes of a CylinderCollection consists primarily of summarry statistics. \
# Statistics may be calculated using dedicated functions, or they may be calculated via the overarching 'statistics' function. \
# (see statistics_calculations for more info)
#

# %%

# Demonstrating several options for working with statistics

# A few statistic specific functions

# Results obtained through the bulk statistics function
