"""Python setup.py for canopyHydrodynamics package"""
from __future__ import annotations

import os

from setuptools import find_packages, setup


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("canhydro", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="canhydro",
    version=read("canhydro", "VERSION"),
    description="Utils for modeling canopy hydrodynamics with QSM's",
    url="https://github.com/wischmcj/canopyHydrodynamics/",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["testopolis", ".github"]),
    install_requires=read_requirements("requirements.txt"),
    entry_points={"console_scripts": ["cylinders = cylinders.__main__:main"]},
)
