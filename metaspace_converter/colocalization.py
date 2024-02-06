from typing import Optional, Tuple

import numpy as np
from anndata import AnnData
from scipy import ndimage

from metaspace_converter.anndata_to_array import anndata_to_image_array
from metaspace_converter.constants import COLOCALIZATION


def coloc_ml_preprocessing(
    adata: AnnData,
    layer: Optional[str] = "coloc_ml_preprocessing",
    median_filter_size: Tuple[int, int] = (3, 3),
    quantile_threshold: float = 0.5,
):
    """
    Preprocessing for colocalization analysis according to the colocML publication
    (https://doi.org/10.1093/bioinformatics/btaa085).

    In the publication the authors evaluated colocalization metrics and preprocessing approaches.
    They found the best performance for
    1) median filtering of ion images with a (3, 3) kernel size and
    2) quantile thresholding ad 50%, meaning all pixels with intensities below the 50%
    quantile set to 0.

    This function performs the same preprocessing steps.
    Recommended for call before running the ``colocalization`` function.

    Args:
        adata: An AnnData object.
        layer: Key for the adata.layers dict in which the processed data will be saved.
            Default value is ``colocml_preprocessing``.
            If None, ``adata.X`` will be overwritten with the processed data.
        median_filter_size: 2-dimensional filter size for the median filtering that
            will be performed per ion image.
        quantile_threshold: Float between 0 and 1. The function computes the quantile value per
            ion image and all pixels below the quantile threshold will be set to 0.

    Returns:
        None. The processed data is saved in ``layer``. If layer is set to None, ``adata.X`` will
        be overwritten

    Raises:
        ValueError: If no annotations are available in ``adata.X``.

    """

    # Extract image array from anndata:
    imarray = anndata_to_image_array(adata=adata)

    if len(imarray) == 0:
        raise ValueError("AnnData contains no annotations. ColocML preprocessing cannot be called.")

    # Median filtering
    imarray = ndimage.median_filter(imarray, size=(1, median_filter_size[0], median_filter_size[1]))

    # reshaping
    imarray = imarray.reshape((imarray.shape[0], -1))

    # Quantile thresholding
    mask = imarray < np.percentile(imarray, q=quantile_threshold * 100, axis=1)[:, np.newaxis]
    imarray[mask] = 0

    # Inserting of processed data into adata object
    if layer is None:
        adata.X = imarray.transpose()
    else:
        adata.layers[layer] = imarray.transpose()


def colocalization(adata: AnnData, layer: Optional[str] = "coloc_ml_preprocessing"):
    """
    Colocalization of ion images using the cosine similarity metric.

    In combination with the ``colocML_preprocessing`` function, this metric performed best in the
    colocML publication (https://doi.org/10.1093/bioinformatics/btaa085).

    It is recommended to call the the ``coloc_ml_preprocessing`` function beforehand.

    Args:
        adata: An AnnData object.
        layer: Key for ``adata.layer`` from which the ionimage_data for preprocessing taken.
            If ``None``, ``adata.X`` is used. ``coloc_ml_preprocessing`` will save the preprocessed
            data per default in ``adata.layer['coloc_ml_preprocessing']``.

    Returns:
        None. The processed data is saved in ``adata.varp['colocalization']``.

    Raises:
        ValueError: If layer is not found in adata.layers.
    """

    # Select data
    if layer is None:
        data = np.array(adata.X).transpose()
    elif layer in adata.layers.keys():
        data = np.array(adata.layers[layer]).transpose()
    else:
        raise ValueError(f"Layer `{layer}` not found in adata.layers.")

    # Compute colocalization
    coloc = _pairwise_cosine_similarity(data)

    # Save colocalization
    adata.varp[COLOCALIZATION] = coloc


def _pairwise_cosine_similarity(data: np.ndarray) -> np.ndarray:

    # Divide image vectors by euclidean norm
    norm_data = data / np.linalg.norm(data, axis=1, keepdims=True)

    # Compute pairwise cosine similarity
    similarity_matrix = np.dot(norm_data, norm_data.T)

    return similarity_matrix
