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

from typing import List

import pytest
from _pytest.nodes import Item

from canhydro.Cylinder import Cylinder
from canhydro.CylinderCollection import CylinderCollection
from canhydro.Forester import Forester
from canhydro.global_vars import test_input_dir


def pytest_collection_modifyitems(items: list[Item]):
    for item in items:
        if "spark" in item.nodeid:
            item.add_marker(pytest.mark.spark)
        elif "_int_" in item.nodeid:
            item.add_marker(pytest.mark.integration)


@pytest.fixture
def basic_cylinder():
    cyl = Cylinder()
    return Cylinder()
    pass


@pytest.fixture
def basic_forest():
    forest = Forester("1_TenCyls.csv")
    forest.get_file_names(dir=test_input_dir)
    forest.qsm_from_file_names()
    return forest
