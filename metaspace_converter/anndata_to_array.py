from typing import Optional
import numpy as np
from anndata import AnnData

from metaspace_converter.constants import COL, METASPACE_KEY, X, Y
from metaspace_converter.to_anndata import all_image_pixel_coordinates


def _check_pixel_coordinates(adata: AnnData) -> bool:
    sorted_obs = adata.obs.sort_values([COL.ion_image_pixel_y, COL.ion_image_pixel_x])

    pixel_list = sorted_obs[[COL.ion_image_pixel_y, COL.ion_image_pixel_x]].values

    img_size = adata.uns[METASPACE_KEY]["image_size"]
    required_pixels = all_image_pixel_coordinates((img_size[Y], img_size[X]))

    return np.all(np.equal(pixel_list, required_pixels))


def anndata_to_image_array(adata: AnnData, layer: Optional[str]=None) -> np.ndarray:
    """
    Extracts an array of ion images from an AnnData object
    (that has been generated through the ``metaspace_to_anndata`` function).

    Args:
        adata: An AnnData object.
        layer: ``AnnData.layer`` that should be extracted to an image array. 
            Default is None, which means that ``adata.X`` will be used.

    Returns:
        A three-dimensional Numpy array in the following shape

            * Dimension 0: Number of ion images in the order of ``adata.var_names``
            * Dimension 1: Image height ``adata.uns["metaspace"]["image_size"]["y"]``
            * Dimension 2: Image width ``adata.uns["metaspace"]["image_size"]["x"]``

    Raises:
        ValueError: If the AnnData object has been modified.
            E.g. Pixel have been removed/added and the number of pixels
            and their coordinates do not match the original image dimensions.
    """
    if layer is None:
        pixel_array = np.array(adata.X).transpose().copy()
    elif layer in adata.layers.keys():
        pixel_array = np.array(adata.layers[layer]).transpose().copy()
    else:
        raise ValueError(f"Layer `{layer}` not found in adata.layers.")
    
    img_size = adata.uns[METASPACE_KEY]["image_size"]

    # Check if image dimensions are okay
    if img_size[X] * img_size[Y] != pixel_array.shape[1]:
        raise ValueError("Number of observations does not match the original image dimensions")

    # Check if all pixels are available
    if not _check_pixel_coordinates(adata):
        raise ValueError("Not all pixels for ion images are available")

    # Sort indices, in case of modified order of pixels (obs)
    image_sorting = adata.obs.sort_values(
        [COL.ion_image_pixel_y, COL.ion_image_pixel_x]
    ).index.values.astype(int)
    pixel_array = pixel_array[:, image_sorting]

    image_array = pixel_array.reshape(
        (
            pixel_array.shape[0],
            adata.obs[COL.ion_image_pixel_y].max() + 1,
            adata.obs[COL.ion_image_pixel_x].max() + 1,
        )
    )

    return image_array
