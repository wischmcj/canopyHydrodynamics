
This page provides an overview of the objects used by canoPyHydro.

# Cylinders
The Cylinder class is used to represent the 3-D cylinders that make up a QSM. The most important function of Cylinder objects is their ability to return data regarding the projections onto planes. Cylider objects utilize our custom 'geometry' module to calculate their projections onto the XY, XZ and YZ planes.

# Cylinder Collections
Cylinder Collections are just as they sound and, at the most basic level, a Cylinder Collection is defined as a list of 1 or more Cylinder objects. Cylinder Collections almost always represent QSM's (or parts of a QSM), and are meant to help users explore their QSMs. Below, we demonstrate how one might initialize a cylinder collection using cylinder data (e.g. QSM data) stored in a CSV file.

# Foresters
Forester objects allow users to conveniently create and manage Cylinder Collections. In particular, Foresters are useful for reading in and processing QSM files. When a Forester object is created, available file names are read from the default directory, './data/input/' or, if specified, from the directory specified by the user. This list of available files can be accessed through the Forester.file_names attribute.
