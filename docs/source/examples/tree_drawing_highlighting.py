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

# %% [markdown]
# # Drawing and Highlighting

# %%
# For the purposes of this tutorial, we will turn off logging 
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

# %% [markdown]
# Let's first discuss our low code approach to figure creation \
# Examine the output produced by the following two cells. Note that their outputs are the same, but the former uses quite a bit more code. 

# %%

from canopyhydro.CylinderCollection import CylinderCollection
# Initializing a CylinderCollection object
myCollection = CylinderCollection()

# Converting a specified file to a CylinderCollection object 
myCollection.from_csv('example_tree.csv')

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.project_cylinders('XZ')
myCollection.draw('XZ')

# Requesting an plot of the tree projected onto the YZ plane ('side' view)
myCollection.project_cylinders('YZ')
myCollection.draw('YZ')

# Requesting an plot of the tree projected onto the XY plane (Birds eye view)
myCollection.project_cylinders('XY')
myCollection.draw('XY')

# %%
# Below We can see the minimal code needed to plot a tree
from canopyhydro.CylinderCollection import CylinderCollection
from matplotlib import pyplot as plt    

# Requesting an plot of the tree projected onto the XZ plane ('fromt' view)
myCollection.draw('XZ')

# Requesting an plot of the tree projected onto the YZ plane ('side' view)
myCollection.draw('YZ')

# Requesting an plot of the tree projected onto the XY plane (Birds eye view)
myCollection.draw('XY')




# %% [markdown]
# The latter, more sucinct option is applocable to not just cylinder projections, but
# for all plotted shapes this includes watershed boundaries (e.g. convex hull), drip points as well as their 'XY,' 'XZ' and 'ZY' projections
# This is to say, if a requested object is not available when the draw function is called, then the object will be calculated with the default attributes and drawn /

# %%

# %%
# Highlighting the stem flow component
from canopyhydro.CylinderCollection import CylinderCollection
from matplotlib import pyplot as plt    

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
