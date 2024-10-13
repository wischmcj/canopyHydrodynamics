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

from __future__ import annotations

from src.canopyhydro.CylinderCollection import CylinderCollection

# %% [markdown]
# # Use Case Example #1
# %% [markdown]
# ## Setup
#
# In this example, we will show case
#

# %% [markdown]
#

# %%


# Initializing a CylinderCollection object
myCollection = CylinderCollection()

# Converting a specified file to a CylinderCollection object
myCollection.from_csv('Secrest10-08_000000.csv')

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.project_cylinders('XZ')
myCollection.draw('XZ')

# Requesting an plot of the tree projected onto the YZ plane ('side' view)
myCollection.project_cylinders('YZ')
myCollection.draw('YZ')

# Requesting an plot of the tree projected onto the XY plane (Birds eye view)
myCollection.project_cylinders('XY')
myCollection.draw('XY')

# %% [markdown]
#

from matplotlib import pyplot as plt

# %%
# Below We can see the minimal code needed to plot a tree
from src.canopyhydro.CylinderCollection import CylinderCollection

# Initializing a CylinderCollection object
myCollection = CylinderCollection()

# Converting a specified file to a CylinderCollection object
myCollection.from_csv('example_tree.csv')

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.draw('XZ')

# Requesting an plot of the tree projected onto the YZ plane ('side' view)
myCollection.draw('YZ')

# Requesting an plot of the tree projected onto the XY plane (Birds eye view)
myCollection.draw('XY')


print('Finished basic tree plot')


from matplotlib import pyplot as plt

# %%
# Highlighting the stem flow component
from src.canopyhydro.CylinderCollection import CylinderCollection

# Initializing a CylinderCollection object
# myCollection = CylinderCollection()

# Converting a specified file to a CylinderCollection object
# myCollection.from_csv('example_tree.csv')

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
# myCollection.project_cylinders('XY')


# myCollection.initialize_digraph_from()
# myCollection.find_flow_components()
# myCollection.calculate_flows()


myCollection.draw('XY', highlight_lambda=lambda:is_stem, save = True, file_name_ext="docs_ex")
myCollection.draw('XZ', highlight_lambda=lambda:is_stem, save = True, file_name_ext="docs_ex")


# %%
# Drawing various branch orders
myCollection.draw('XZ', filter_lambda=lambda:branch_order<=1, highlight_lambda=lambda:branch_order==1, save = True, file_name_ext="docs_ex")
myCollection.draw('XZ', filter_lambda=lambda:branch_order<=2, highlight_lambda=lambda:branch_order==2, save = True, file_name_ext="docs_ex")
myCollection.draw('XZ', filter_lambda=lambda:branch_order<=3, highlight_lambda=lambda:branch_order==3, save = True, file_name_ext="docs_ex")
myCollection.draw('XZ', filter_lambda=lambda:branch_order<=4, highlight_lambda=lambda:branch_order==4, save = True, file_name_ext="docs_ex")
