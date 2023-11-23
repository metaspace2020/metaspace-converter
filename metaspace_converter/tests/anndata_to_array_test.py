from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

if TYPE_CHECKING:
    import _pytest.fixtures

from metaspace_converter import anndata_to_image_array
from metaspace_converter.constants import COL, METASPACE_KEY, X, Y
from metaspace_converter.to_anndata import all_image_pixel_coordinates


@pytest.fixture
def adata_with_coordinates_and_image_size(
    request: "_pytest.fixtures.SubRequest",
) -> tuple[AnnData, np.ndarray]:
    """
    Create an AnnData object filled with parametrized shape

    Args:
        request: Request passed by Pytest, can be parametrized with a dictionary containing
            ``num_ions``, ``height``, ``width``

    Returns:
        A tuple of an AnnData object and the expected array of stacked ion images
    """
    num_ions = request.param.get("num_ions", 0)
    height = request.param.get("height", 0)
    width = request.param.get("width", 0)
    shape = (num_ions, height, width)
    # Create a stack of ion images where every pixel has a different value
    ion_images_stack = np.arange(np.prod(shape)).reshape(shape)
    coordinates = all_image_pixel_coordinates(shape[1:])
    # Create an AnnData with pixel coordinates
    obs = pd.DataFrame(
        {COL.ion_image_pixel_y: coordinates[:, 0], COL.ion_image_pixel_x: coordinates[:, 1]}
    )
    adata = AnnData(
        X=ion_images_stack.reshape((height * width, num_ions)),
        obs=obs,
        uns={METASPACE_KEY: {"image_size": {X: width, Y: height}}},
    )
    expected = ion_images_stack
    return adata, expected


@pytest.mark.parametrize(
    "adata_with_coordinates_and_image_size",
    [
        # Non-square ion images
        dict(num_ions=4, height=2, width=3),
        # Edge case: No annotations found
        dict(num_ions=0, height=2, width=3),
    ],
    indirect=["adata_with_coordinates_and_image_size"],
)
def test_anndata_to_image_array(adata_with_coordinates_and_image_size: AnnData):
    adata, expected = adata_with_coordinates_and_image_size

    actual = anndata_to_image_array(adata)

    assert actual.shape == (
        adata.shape[1],
        adata.obs[COL.ion_image_pixel_y].max() + 1,
        adata.obs[COL.ion_image_pixel_x].max() + 1,
    )

    assert actual.shape == (
        adata.shape[1],
        adata.uns[METASPACE_KEY]["image_size"][Y],
        adata.uns[METASPACE_KEY]["image_size"][X],
    )
