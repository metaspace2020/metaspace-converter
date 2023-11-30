.. METASPACE Converter documentation master file, created by
   sphinx-quickstart on Mon Jan 30 10:10:43 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to METASPACE Converter's documentation!
===============================================

|Tests badge|
|Docs badge|
|PyPI badge|

Python package to download and convert datasets from the `METASPACE`_
knowledge base to common formats for single-cell and spatial omics analysis.
Datasets can be directly downloaded to
`AnnData`_ and `SpatialData`_ objects.

`AnnData`_ is the underlying data format
of many packages of the `scverse`_ such as
`ScanPy`_ for single-cell data analysis and
`SquidPy`_ for spatial omics analysis.

Another supported format that is part of the `scverse`_
is `SpatialData`_ for storing, aligning, and processing spatial omics data. 
This enables users to easily align and integrate METASPACE datasets to other spatial omics modalities.

The METASPACE-converter package uses the `METASPACE python client`_
download datasets from METASPACE.
If you also need to upload or modify datasets on METASPACE, please check the `Python client documentation`_.

If you encounter bugs or have suggestions for new features, please open an issue on `GitHub`_.


Installation
------------
You can install the package via pip:


.. code-block::
   
   pip install metaspace-converter 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api
   examples

License
-------

Unless specified otherwise source code in file headers or LICENSE files present in subdirectories, all files of this package are licensed under the `Apache 2.0 license`_.



.. _METASPACE: https://metaspace2020.eu/
.. _AnnData: https://anndata.readthedocs.io/en/stable/
.. _ScanPy: https://scanpy.readthedocs.io/en/stable/
.. _SquidPy: https://squidpy.readthedocs.io/en/stable/
.. _SpatialData: https://spatialdata.scverse.org/en/latest/
.. _scverse: https://doi.org/10.1038/s41587-023-01733-8
.. _GitHub: https://github.com/metaspace2020/metaspace-converter
.. _METASPACE python client: https://github.com/metaspace2020/metaspace/tree/master/metaspace/python-client
.. _Python client documentation: https://metaspace2020.readthedocs.io/en/latest/index.html
.. _Apache 2.0 license: https://github.com/metaspace2020/metaspace/blob/master/LICENSE

.. |Tests badge| image:: https://img.shields.io/github/actions/workflow/status/metaspace2020/metaspace-converter/tests.yml?branch=master&label=tests
    :target: https://github.com/metaspace2020/metaspace-converter/actions/workflows/tests.yml
.. |Docs badge| image:: https://img.shields.io/github/actions/workflow/status/metaspace2020/metaspace-converter/docs.yml?label=documentation
    :target: https://metaspace2020.github.io/metaspace-converter/
.. |PyPI badge| image:: https://img.shields.io/pypi/v/metaspace-converter
    :target: https://pypi.org/project/metaspace-converter/
