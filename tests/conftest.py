#   ---------------------------------------------------------------------------------
#   Copyright (c) Microsoft Corporation. All rights reserved.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   ---------------------------------------------------------------------------------
"""
This is a configuration file for pytest containing customizations and fixtures.

In VSCode, Code Coverage is recorded in config.xml. Delete this file to reset reporting.
"""

# https://docs.pytest.org/en/7.1.x/example/markers.html
# This page outlines the marks available to us
from __future__ import annotations

import pytest
from _pytest.nodes import Item

from canhydro.Forester import Forester
from canhydro.global_vars import test_input_dir


def pytest_collection_modifyitems(items: list[Item]):
    for item in items:
        if "spark" in item.nodeid:
            item.add_marker(pytest.mark.spark)
        elif "_int_" in item.nodeid:
            item.add_marker(pytest.mark.integration)


# @lru_cache(maxsize=256)
@pytest.fixture
def basic_forest():
    forest = Forester()
    forest.get_file_names(dir=test_input_dir)
    # forest.qsm_from_file_names()
    return forest


# @lru_cache(maxsize=256)
@pytest.fixture
def flexible_collection(basic_forest, request):
    basic_forest.qsm_from_file_names(file_name=request.param)
    flexible_collection = basic_forest.cylinder_collections[0]
    flexible_collection.project_cylinders("XZ")
    flexible_collection.project_cylinders("XY")
    flexible_collection.initialize_graph()
    return flexible_collection


@pytest.fixture
def ez_projection():
    forest = Forester()
    forest.get_file_names(dir=test_input_dir)
    forest.qsm_from_file_names(file_name="2_EZ_projection.csv")
    collection = forest.cylinder_collections[0]
    return collection


@pytest.fixture
def happy_path_projection():
    forest = Forester()
    forest.get_file_names(dir=test_input_dir)
    forest.qsm_from_file_names(file_name="3_HappyPathWTrunk.csv")
    collection = forest.cylinder_collections[0]
    return collection


@pytest.fixture
def small_tree():
    forest = Forester()
    forest.get_file_names(dir=test_input_dir)
    forest.qsm_from_file_names(file_name="5_SmallTree.csv")
    collection = forest.cylinder_collections[0]
    return collection


@pytest.fixture
def large_collection():
    forest = Forester()
    forest.get_file_names(dir=test_input_dir)
    forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
    collection = forest.cylinder_collections[0]
    # collection.project_cylinders("XZ")
    return collection
