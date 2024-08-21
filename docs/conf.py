# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from __future__ import annotations

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
project = "canoPyHydro"
copyright = '2024, [{name = "Collin Wischmeyer", email = "cjwischmeyer@gmail.com"}, {name = "Travis Swanson", email = "travis.swanson@gmail.com"}, {name = "John Van Stan", email = "j.vanstan@csuohio.edu"},]'
author = '[{name = "Collin Wischmeyer", email = "cjwischmeyer@gmail.com"}, {name = "Travis Swanson", email = "travis.swanson@gmail.com"}, {name = "John Van Stan", email = "j.vanstan@csuohio.edu"},]'
release = "0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = [
    "myst_parser",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]
