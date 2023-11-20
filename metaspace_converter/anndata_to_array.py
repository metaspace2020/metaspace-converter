import numpy as np
from anndata import AnnData


def _check_pixel_list(adata: AnnData) -> bool:
    
    sorted_obs = adata.obs.sort_values(['ion_image_pixel_y', 'ion_image_pixel_x'])
    
    pixel_list = [(y, x) for y, x in zip(sorted_obs['ion_image_pixel_y'], sorted_obs['ion_image_pixel_x'])]
    
    img_size = adata.uns['metaspace']['image_size']
    required_pixels = [(y, x) for y in range(img_size['y']) for x in range(img_size['x'])]

    return pixel_list == required_pixels


def anndata_to_image_array(adata: AnnData) -> np.ndarray:
    """
    Extracts an array of ion images from an AnnData object 
    (that has ben generated through the `metaspace_to_anndata` function). 

    Args:
        adata: An AnnData object.

    Returns:
        A three-dimensional numpy array in the following shape
        - Dimension 0: Number of ion images in the order of adata.var.shape[0]
        - Dimension 1: Image height (adata.uns['metaspace']['image_size']['y']
        - Dimension 2: Image width (adata.uns['metaspace']['image_size']['x']

    Raises:
        ValueError: If the AnnData object has been modified. E.g. Pixel have been removed/added and the number of pixels 
        and their coordinates do not match the original image dimensions.
    """
    
    pixel_array = adata.X.transpose().copy()
    img_size = adata.uns['metaspace']['image_size']
    
    # Check if image dimensions are okay
    if img_size['x']*img_size['y'] != pixel_array.shape[1]:
        raise ValueError('Number of observations does not match the original image dimensions')
    
    # Check if all pixels are available
    if not _check_pixel_list(adata):
        raise ValueError('Not all pixels for ion images are available')
    
    # Sort incides, in case of modified order of pixels (obs)
    image_sorting = adata.obs.sort_values(['ion_image_pixel_y', 'ion_image_pixel_x']).index.values.astype(int)
    pixel_array = pixel_array[:, image_sorting]
    
    image_array = pixel_array.reshape((pixel_array.shape[0], 
                                       adata.obs['ion_image_pixel_y'].max()+1,
                                       adata.obs['ion_image_pixel_x'].max()+1
                                      ))
        
    return image_array
