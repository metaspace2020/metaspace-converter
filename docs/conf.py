# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "METASPACE Converter"
copyright = "2023; Tim Daniel Rose, Andreas Eisenbarth, Theodore Alexandrov"
author = "Tim Daniel Rose, Andreas Eisenbarth"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # Include API documentation from docstrings
    "sphinx.ext.autodoc",
    # Executes and tests example code to ensure it is (still) working
    "sphinx.ext.doctest",
    # Links symbols in docstrings to API documentation of external projects
    "sphinx.ext.intersphinx",
    # Support for Google style docstrings
    "sphinx.ext.napoleon",
    # Adds a source code page to each API element
    "sphinx.ext.viewcode",
    # Uses type hints from signature instead of types in docstrings
    "sphinx_autodoc_typehints",
    # Links symbols within code blocks to API documentation
    "sphinx_codeautolink",
    # Adds a copy button to code blocks
    "sphinx_copybutton",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Doctest allows to execute example code with directives "" using "make doctest".
# Compared to the Python standard library module, this extension detects code blocks
# without requiring ">>> " as prefix.
# Remove doctest flags within the examples ("#doctest: +SKIP"), not to distract readers.
trim_doctest_flags = True
# Setup common imports, so they are available in the examples.
doctest_global_setup = """
try:
    import numpy as np
    import pandas as pd
    import matplotlib
    # Ensure doctest does not attempt to display plots in blocking GUI windows.
    matplotlib.use("agg")
except ImportError:
    np = None
    pd = None
"""

# This creates links of symbols to external documentation.
intersphinx_mapping = {
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "metaspace": ("https://metaspace2020.readthedocs.io/en/latest/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
    "squidpy": ("https://squidpy.readthedocs.io/en/stable/", None),
    "spatialdata": ("https://spatialdata.scverse.org/en/latest/", None),
    "spatialdata-plot": ("https://spatialdata.scverse.org/projects/plot/en/latest/", None),
}
# Needed to help MyST distinguish between URL schemes and references with a colon (?)
myst_url_schemes = ["http", "https"]
