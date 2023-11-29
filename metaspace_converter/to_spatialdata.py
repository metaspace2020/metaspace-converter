import warnings
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from metaspace import SMInstance
from metaspace.sm_annotation_utils import SMDataset
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel, TableModel
from spatialdata.transformations import Affine, Scale, Sequence as SequenceTransform, Translation

from metaspace_converter.constants import (
    COL,
    COORD_SYS_GLOBAL,
    INSTANCE_KEY,
    MICROMETER,
    OPTICAL_IMAGE_KEY,
    POINTS_KEY,
    REGION_KEY,
    XY,
    YX,
    YXC,
    X,
    Y,
)
from metaspace_converter.to_anndata import metaspace_to_anndata
from metaspace_converter.utils import has_optical_image


def metaspace_to_spatialdata(
    dataset: Optional[SMDataset] = None,
    dataset_id: Optional[str] = None,
    database: Optional[tuple[str, str]] = None,
    fdr: float = 0.1,
    use_tic: bool = False,
    metadata_as_obs: bool = False,
    add_optical_image: bool = True,
    optical_name_added: str = OPTICAL_IMAGE_KEY,
    add_points: bool = True,
    points_name_added: str = POINTS_KEY,
    sm: Optional[SMInstance] = None,
    **annotation_filter,
) -> SpatialData:
    """
    Download a METASPACE dataset as a :py:class:`SpatialData` object.

    See: https://metaspace2020.eu/about

    Args:
        dataset: A METASPACE dataset instance. If not provided, the ``dataset_id`` must be given.
        dataset_id: The unique ID of a dataset on METASPACE, e.g. "2021-09-03_11h43m13s"
        database: A single METASPACE database given as a tuple of name and version. Usually it is
            displayed on METASPACE as "HMDB – v4" which corresponds to ``("HMDB", "v4")``.
        fdr: Returns only annotations for which the false discovery rate is less or equal to this
            limit.
        use_tic: When True, the output values will be scaled by the total ion count per pixel and
            will be in 0.0 to 1.0 range.
        metadata_as_obs: Whether to store metadata in the ``obs`` dataframe instead of ``uns``. For
            a single METASPACE dataset, metadata is the same for all pixels, so it would be
            duplicated for all ``obs``. When combining multiple datasets, it would be preserved in
            ``obs`` but not in ``uns``.
        add_optical_image: Whether to also add the optical image (if one exists) to SpatialData
            images. If none exists, it is not added, and no error raised.
        optical_name_added: Name of the element where to store the image in the SpatialData object
        add_points: Whether to also add ion image pixel coordinates as SpatialData points. This
            allows to spatially visualize the ion image values. If False, only ion intensities are
            added to the table and coordinates are added to ``obs``.
        points_name_added: Name of the element where to store the points in the SpatialData object
        sm: Optionally a cached SMInstance
        annotation_filter: Additional keyword arguments passed to the METASPACE API.

    Returns:
        A :py:class:`SpatialData` object with
            * ion intensities and metadata: ``.table``
            * optical image: ``.images["optical_image"]`` (or ``optical_name_added``)
            * sampling coordinates: ``.points["maldi_points"]`` (or ``points_name_added``) in a
              coordinate system relative to the top-left image corner and physical scale
              (micrometers)

    Raises:
        ValueError: If something is wrong with the input data or parameters
        metaspace.sm_annotation_utils.DatasetNotFound: If the dataset ID is not found
        ConnectionError: For any other network connectivity problems
    """
    if dataset is None:
        if sm is None:
            sm = SMInstance()
        dataset = sm.dataset(id=dataset_id)

    # Create table
    adata = metaspace_to_anndata(
        dataset=dataset,
        database=database,
        fdr=fdr,
        use_tic=use_tic,
        metadata_as_obs=metadata_as_obs,
        add_optical_image=False,
        sm=sm,
        **annotation_filter,
    )
    table = _create_spatialdata_table(adata)
    sdata = SpatialData(table=table)

    # Create image
    if has_optical_image(dataset) and add_optical_image:
        optical_image = optical_image_to_spatial_image(dataset=dataset)
        sdata.images[optical_name_added] = optical_image

    # Create points
    if add_points:
        points = _create_points(dataset, adata)
        sdata.points[points_name_added] = points

    return sdata


def _is_matrix_homogeneous(transform_matrix_2d: np.ndarray) -> bool:
    # In an n×n matrix, the last row (except last element) should be zeros.
    actual = transform_matrix_2d[-1, :-1]
    expected = np.zeros(transform_matrix_2d.shape[1] - 1)
    return np.array_equal(actual, expected)


def optical_image_to_spatial_image(dataset: SMDataset) -> Optional[SpatialImage]:
    """
    Download the optical image of a METASPACE dataset to a SpatialImage object.

    This function requires that the dataset has an optical image. Otherwise, an error is raised.

    Args:
        dataset: A METASPACE dataset

    Returns:
        A SpatialImage object with the image's transformation to world coordinates
    """
    optical_images = dataset.optical_images()
    image_yx_rgb = optical_images[0]
    # Transformation from optical image to ion image coordinate system (inverse of ``_transforms``)
    matrix_xy = optical_images._itransforms[0]
    if _is_matrix_homogeneous(matrix_xy):
        affine_matrix_xy = matrix_xy
    else:
        warnings.warn(
            "The METASPACE optical image has a projective transformation, but "
            "SpatialData only supports affine transformations. "
            "Using an approximated affine."
        )
        affine_matrix_xy = matrix_xy.copy()
        affine_matrix_xy[-1, :-1] = 0.0
    optical_to_ion_image = Affine(matrix=affine_matrix_xy, input_axes=XY, output_axes=XY)
    # Transformation from ion image to world coordinate system
    ion_image_to_global = get_ion_image_to_physical_transformation(dataset)
    pixel_corner_to_center = Translation(translation=[-0.5, -0.5], axes=YX)
    optical_to_global = SequenceTransform(
        [optical_to_ion_image, pixel_corner_to_center, ion_image_to_global]
    ).to_affine(input_axes=YX, output_axes=YX)
    return Image2DModel.parse(
        image_yx_rgb,
        transformations={COORD_SYS_GLOBAL: optical_to_global},
        dims=YXC,
        scale_factors=(2, 2, 2),
        axis_units={Y: MICROMETER, X: MICROMETER},
        c_coords=("r", "g", "b"),
        # rgb=True,
    )


def _create_spatialdata_table(adata: AnnData) -> AnnData:
    # Add the index as a separate column for matching rows to regions in labels or points.
    adata.obs[INSTANCE_KEY] = adata.obs.index.astype(int)
    # Avoid same name for index and a column name, which is ambiguous and causes warning or error.
    adata.obs.index.name = "index"
    adata.obs[REGION_KEY] = pd.Series(
        [POINTS_KEY] * adata.n_obs, index=adata.obs.index, dtype="category"
    )
    return TableModel.parse(
        adata,
        region=[POINTS_KEY],  # element in SpatialData to be used
        region_key=REGION_KEY,  # table column containing the regions
        instance_key=INSTANCE_KEY,  # table column to be matched to points index
    )


def _create_points(dataset: SMDataset, adata: AnnData) -> pd.DataFrame:
    # FIXME: PointsModel does not yet allow to specify units, we should set micrometers as in image.
    # FIXME: Workaround for that napari-spatialdata does not allow to visualize values from ``X``
    #  array, only from ``obs``. For now, we copy X to ``obs`` to enable visualization.
    # points = PointsModel.parse(
    #     adata.obs[[X, Y]].reset_index(),
    #     transformations={COORD_SYS_GLOBAL: get_ion_image_to_physical_transformation(dataset)},
    # )
    return PointsModel.parse(
        pd.concat(
            [
                adata.obs[[COL.ion_image_pixel_x, COL.ion_image_pixel_y]].reset_index(),
                pd.DataFrame(adata.X, columns=adata.var.index),
            ],
            axis=1,
        ),
        coordinates={"x": COL.ion_image_pixel_x, "y": COL.ion_image_pixel_y},
        transformations={COORD_SYS_GLOBAL: get_ion_image_to_physical_transformation(dataset)},
    )


def get_ion_image_to_physical_transformation(dataset: SMDataset) -> Scale:
    """
    Get the transformation from the ion image (pixels) to the world coordinate system (micrometers).

    Args:
        dataset: A METASPACE dataset

    Returns:
        A scale SpatialData transformation
    """
    ion_image_pixel_size = dataset.metadata["MS_Analysis"]["Pixel_Size"]
    scale_factors_xy = [ion_image_pixel_size["Xaxis"], ion_image_pixel_size["Yaxis"]]
    return Scale(scale=scale_factors_xy, axes=XY)
