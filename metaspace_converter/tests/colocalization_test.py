from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

if TYPE_CHECKING:
    import _pytest.fixtures

from metaspace_converter import anndata_to_image_array
from metaspace_converter.colocalization import colocalization, coloc_ml_preprocessing
from metaspace_converter.constants import COL, COLOCALIZATION, METASPACE_KEY, X, Y
from metaspace_converter.to_anndata import all_image_pixel_coordinates


@pytest.fixture
def adata_dummy(request: "_pytest.fixtures.SubRequest") -> AnnData:
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
    return adata


@pytest.mark.parametrize(
    "adata_dummy",
    [
        dict(num_ions=5, height=10, width=10),
        dict(num_ions=0, height=10, width=10),
        dict(num_ions=5, height=3, width=4),
    ],
    indirect=["adata_dummy"],
)
def test_colocml_preprocessing(adata_dummy):

    COLOCML_LAYER = "coloc_ml_preprocessing"

    adata = adata_dummy

    if adata.X.shape[1] == 0:
        with pytest.raises(ValueError) as e_info:
            coloc_ml_preprocessing(
                adata, median_filter_size=(3, 3), quantile_threshold=0, layer=COLOCML_LAYER
            )
    else:
        

        # Test median filtering
        coloc_ml_preprocessing(
            adata, median_filter_size=(3, 3), quantile_threshold=0, layer=COLOCML_LAYER
        )

        actual = anndata_to_image_array(adata)

        # Layer exists
        assert COLOCML_LAYER in adata.layers.keys()

        # Layer sizes match
        assert adata.X.shape == adata.layers[COLOCML_LAYER].shape

        # median filter
        # Check that preprocessing was done correctly.
        # Probe a single pixel of the resulting array and compute the expected median filter.
        x = 1
        y = 1
        width = adata.uns[METASPACE_KEY]["image_size"][X]
        expected_pixel_features = np.median(actual[:, y - 1:y + 2, x - 1:x + 2], axis=[1, 2])
        actual_pixel_features = adata.layers[COLOCML_LAYER][y * width + x, :]
        np.testing.assert_allclose(actual_pixel_features, expected_pixel_features)

        # expected_preprocessing = np.median(actual[0][:3, :3])
        # observed = adata.layers[COLOCML_LAYER][adata.uns[METASPACE_KEY]["image_size"][X] + 1, 0]
        # assert observed == expected_preprocessing

        # expected_preprocessing = np.median(actual[1][:3, :3])
        # observed = adata.layers[COLOCML_LAYER][adata.uns[METASPACE_KEY]["image_size"][X] + 1, 1]
        # assert observed == expected_preprocessing

        # Quantile thresholding

        coloc_ml_preprocessing(
            adata, median_filter_size=(3, 3), quantile_threshold=0.5, layer=COLOCML_LAYER
        )

        assert all(np.sum(adata.layers[COLOCML_LAYER] == 0, axis=0) <= 0.5 * adata.X.shape[0])

       


def test_colocalization():

    arr = np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [1.0, 0.0, 0.0]]).transpose()

    expected = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    adata = AnnData(X=arr)

    colocalization(adata, layer=None)

    assert COLOCALIZATION in adata.varp.keys()

    coloc = adata.varp[COLOCALIZATION]

    assert np.all(coloc == expected)
