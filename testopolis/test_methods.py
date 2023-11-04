"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import pytest

from src import Cylinder
from src import CylinderCollection
from src import Forester

def test_hello(unit_test_mocks: None):
    """
    This is a simple test, which can use a mock to override online functionality.
    unit_test_mocks: Fixture located in conftest.py, implictly imported via pytest.
    """
    hello_test()

def test_split(self):
    s = 'hello world'
    """This is an example of a type error test"""
    pytest.assertEqual(s.split(), ['hello', 'world'])
    # check that s.split fails when the separator is not a string
    with self.assertRaises(TypeError):
        s.split(2)

def test_file_names(self, dir:str):
    forest = Forester( directory=dir )
    print(forest.getFileNames())
    print(forest._file_names)
    self.assertEqual(forest._file_names, ['4_Unhappy_DripPathAdjTrunkWTrunk.csv', '2A_DripPathBegStartMultHappy.csv', '1_HappyPathWTrunk.csv'])
    print('File Names Successfull')

def test_create_cyliders(self):
    forest = Forester( directory=testDir )
    forest.getFileNames()
    
    forest.QSM_FromFileNames()
    return forest
    # self.assertTrue('FOO'.isupper())