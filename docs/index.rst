.. METASPACE Converter documentation master file, created by
   sphinx-quickstart on Mon Jan 30 10:10:43 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to METASPACE Converter's documentation!
===============================================

|Tests badge|
|Docs badge|
|PyPI badge|

Functions to convert `METASPACE`_ datasets to `AnnData`_ objects.

This makes it easy to work with the `ScanPy`_, `SquidPy`_ and `SpatialData`_ ecosystem with METASPACE data.

Installation
------------
You can install the package via pip:

.. testcode::
   
   pip install metaspace-converter 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api
   examples

.. _METASPACE: https://metaspace2020.eu/
.. _AnnData: https://anndata.readthedocs.io/en/stable/
.. _ScanPy: https://scanpy.readthedocs.io/en/stable/
.. _SquidPy: https://squidpy.readthedocs.io/en/stable/
.. _SpatialData: https://spatialdata.scverse.org/en/latest/


.. |Tests badge| image:: https://img.shields.io/github/actions/workflow/status/metaspace2020/metaspace-converter/tests.yml?branch=master&label=tests
    :target: https://github.com/metaspace2020/metaspace-converter/actions/workflows/tests.yml
.. |Docs badge| image:: https://img.shields.io/github/actions/workflow/status/metaspace2020/metaspace-converter/docs.yml?label=documentation
    :target: https://metaspace2020.github.io/metaspace-converter/
.. |PyPI badge| image:: https://img.shields.io/pypi/v/metaspace-converter
    :target: https://pypi.org/project/metaspace-converter/
