
# Example notebooks

Rather than strictly providing definitions and usage examples, as in the `How To` section, several Jupyter notebooks have been provided here to provide an interactive turotial experience.

These notebooks are packaged with the `canoPyHydro` package, and will therefore be cloned along side the rest of the repo when using git. The notebooks can also be downloaded individually from the `code_samples/notebooks` directory:

```bash
git clone https://github.com/wischmcj/canopyHydrodynamics.git
cd canopyHydrodynamics/docs/source/examples/
```

Note that if you're running locally, in order to run the code in the workflows, you will also need to install the canoPyHydro package from pip. Please refer to the [getting started page](https://canopyhydrodynamics.readthedocs.io/en/latest/getting_started.html) for help getting started.

# Feature Descriptions

**1. [Introduction to Common Objects](intro_to_objects.ipynb)**: This notebook provides a detailed walkthrough of the various different objects used in the canoPyHydro package.

**2. [Projecting Cylinders](projecting_cylinders.ipynb)**: This notebook demonstrates how 3D cylinders that make up a QSM model are converted to 2D shapes via projection.

**3. [Concave Hulls and Watersheds](alpha_shapes.md)**: This notebook demonstrates how 3D cylinders that make up a QSM model are converted to 2D shapes via projection.

# Feature Application

**1. [Drawing and Highlighting](tree_drawing_highlighting.ipynb)**: This notebook focuses on non-flow related visualizations of tree canopies, showcasing canopyHydro's abillity to filter out portions of a tree, highlight features of interest and provide various viewing angles of the results.

**2. [Flow Identification and Drawing](flow_identification_drawing.ipynb)**: This notebook shows various code snippets demonstrating how to use the `canoPyHydro` package to model canopy hydrodynamics. This is inclusive of both: finding stem flow and through fall statistics and visualizing these flows in the tree canopy.
