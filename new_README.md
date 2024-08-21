<!-- LOGO -->

<br />
<h1>
<p align="center">
    <img src=".canopyhydro_logo.jpeg" height="390" width="390">
</p>

  CanoPyMap

</h1>
  <p align="center">
    Leveraging remote sensing to characterize canopy watersheds
    <br />
    </p>
</p>
<p align="center">
  <a href="#about-the-project">About The Project</a> •
  <a href="#usage">How To Use</a> •
  <a href="#examples">Examples</a> •
  <a href="#best-practice">Best Practice</a> •
  <a href="#credits">Credits</a> •
  <a href="examples.md">More Examples</a>
</p>

# canoPyHydro

## Introduction

This package provides functionality to assist researchers in investigating water partitioning in tree canopies. A significant portion of this work is centered around the identification and summarization of canopy structural traits. As such the functionality available in this package includes the following:

## Overview

The methods provided in canoPyHydro extract data from these QSMs (i.e. cylinder angle, projected area, ...) and use said data to predict how percipitation falling in the vicinity of the corresponding tree will be paritioned. At present, these paritions are as follows

1. Stem flow - Percipitation intercepted by the tree that flows down to the trunk of the tree
2. Drip flow - Percipitation intercepted by the tree that drips to the ground prior to reaching the trunk of the tree
3. Throughfall - Percipitation that falls in the gaps of the tree canopy, reaching the ground without coming into contact with the tree

Methods in the former category primarily utilize the orientation of the cylinders - in tandem with user defined parameters - to determine the direction in which intercepted water will flow. Once all cylinders have been assigned to a partition, methods in the latter category can then be applied to provide a wide array of statistics describing each partition - most importantly: contributing branch surface area, 2D projected area and canopy coverage area.

<exerpt from johns doc>

For more details of the methodology please refer to [`common_citations.md`](./common_citations.md).

## Documentation

A complete documentation for the package can be found at the `<instert johns doc here>`

## Installation

First, Create  conda environment

```
python3 -m venv venv
source venv/bin/activate
```

The package is available on PyPI and can be installed using pip

```
pip install canoPyHydro
```

For the source code you can clone the git repository locally using

```
git clone https://github.com/wischmcj/canopyHydrodynamics
```

and then install the dependencies locally using pip

```
cd canopyHydrodynamics
python -m pip install .
```

We can test the environment build by running

```
 pip install -r requirements_dev.txt
pytest
```

For a successful build we should have all tests pass. If this is the case then enjoy using canoPyHydro!

## Examples

We provide several examples to illustrate the tasks that canoPyHydro can perform in [`examples`](./examples). These examples are provided as both an annotated Jupyter notebook or a python script (the scripts are further separated for ease), and each example has a detailed description of its content within the README. `example_function` as a [notebook](./examples/notebooks/example_function.ipynb) or [scripts](./examples/scripts/example_function) are the best place to start for an introduction to the methodology, where we apply it to some simple test functions.

## License

canoPyHydro is an open-source software licensed under the MIT License. Check the details in the [`LICENSE`](./LICENSE) file.

## Citations

If you use this package please cite it appropriately using the 'Cite this repository' dropdown in the right sidebar. Moreover, we also provide a bibliography of previous relevant work in [`common_citations.md`](./common_citations.md). This
