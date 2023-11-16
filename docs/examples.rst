Examples
========

Plain AnnData with ScanPy
-------------------------

.. testcode::

   from metaspace_converter import metaspace_to_anndata
   import scanpy as sc

   # Download data and convert to an AnnData object
   adata = metaspace_to_anndata(
       dataset_id="2022-08-05_17h28m56s",
       fdr=0.1,
       database=("BraChemDB", "2018-01"),
   )

   # Visualization with ScanPy
   sc.pl.spatial(
       adata,
       # Choose first ion for visualization
       color=adata.var.index[0],
       img_key=None,
       spot_size=1,
   )

.. image:: ./_static/img/example_img_sc.png
   :alt: Visualization with ScanPy

SquidPy
-------

.. testcode::

   from metaspace_converter import metaspace_to_anndata
   import squidpy as sq

   # Download dataset with optical background image
   adata = metaspace_to_anndata(
       dataset_id="2022-08-05_17h28m56s",
       fdr=0.1,
       database=("BraChemDB", "2018-01"),
       add_optical_image=True,
   )

   sq.pl.spatial_scatter(
       adata, color=adata.var.index[0], shape="square", img=True, size=15, alpha=0.5
   )

.. image:: ./_static/img/example_img_sq.png
   :alt: Visualization with SquidPy

SpatialData
-----------

Here using a reversed colormap which better represents intense values on bright background.

.. testcode::

   from metaspace_converter import metaspace_to_spatialdata
   import spatialdata_plot  # noqa: Not directly used but extends spatialdata

   # Download dataset with optical background image
   sdata = metaspace_to_spatialdata(
       dataset_id="2022-08-05_17h28m56s",
       fdr=0.1,
       database=("BraChemDB", "2018-01"),
   )

   # Workaround: spatialdata-plot currently does not use points transformation
   sdata.points["maldi_points"] = sdata.transform_element_to_coordinate_system(
       sdata.points["maldi_points"], "global"
   )

   (
       sdata.pl.render_images("optical_image")
       .pl.render_points(
           "maldi_points",
           color=sdata.table.var.index[0],
           alpha=1,
           size=2,
           cmap="viridis_r",
       )
       .pl.show(title=sdata.table.var.index[0], coordinate_systems="global")
   )

.. testoutput::
   :hide:

   ...
