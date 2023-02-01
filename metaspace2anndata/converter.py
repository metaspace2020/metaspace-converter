from metaspace.sm_annotation_utils import SMDataset, IsotopeImages
from anndata import AnnData
from typing import Tuple, List, Optional
from .utils import *


def dataset_to_anndata(ds: SMDataset,
                       database: Optional[Tuple[str, str]] = ("HMDB", "v4"),
                       fdr: Optional[float] = 0.1,
                       results: Optional[pd.DataFrame] = None,
                       all_annotation_images: Optional[List[IsotopeImages]] = None) -> AnnData:
    """
    Convert METASPACE datasets to AnnData object.

    Requires only a dataset (SMDataset), database (e.g. `("HMDB", "v4")`), and fdr.
    Output of `SMDataset.results` and `SMDataset.all_ion_images` can be provided to work,
    since the function does not allow to specify all parameters.

    This function does not download the optical image for datasets.

    Metadata is stored in `adata.uns`.

    :param ds: METASPACE dataset object
    :param database: Database to download results. Is not used if `results` and `all_annotation_images` is provided.
    :param fdr: FDR threshold. Is not used if `results` and `all_annotation_images` is provided.
    :param results: Output of `SMDataset.results()` call.
    :param all_annotation_images: Output of `SMDataset.all_ion_images()` call.
    :return: AnnData object including images, metadata, and annotation information.
    """

    if results is None:
        res = ds.results(database=database, fdr=fdr)
    else:
        res = results.copy()

    res = m_results_replace_index(res)

    if all_annotation_images is None:
        ann_images = ds.all_annotation_images(fdr=fdr,
                                              database=database,
                                              only_first_isotope=True)
    else:
        ann_images = all_annotation_images.copy()

    # Sort them according to res
    img_dict = {image.formula + image.adduct: image for image in ann_images}

    im_array = np.array([img_dict[ion] for ion in res.index])

    # Remove additional dimension with placeholder for isotopes
    im_array = reshape_im_array(im_array)

    # Transform to X,Y representation
    xdata = im_array.reshape((im_array.shape[0], -1)).transpose()

    xy_coord, xy_coord_df = get_xy_coordinates(im_array, xdata)

    adata = AnnData(X=xdata, var=res, obs=xy_coord_df, obsm={"spatial": xy_coord})

    # Only important for background image, but Squidpy sometimes complains, if this is not there.
    spatial_key = "spatial"
    library_id = "image"
    adata.uns[spatial_key] = {library_id: {}}
    adata.uns[spatial_key][library_id]["images"] = {}
    # adata.uns[spatial_key][library_id]["images"] = {"hires": oi[0]}
    adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 1, "spot_diameter_fullres": 0.5}

    # Get METASPACE infos
    add_metadata(adata, ds)

    return adata


def dataset_to_anndata_oi(ds: SMDataset,
                          database: Tuple[str, str] = ("HMDB", "v4"),
                          fdr: float = 0.1,
                          results: pd.DataFrame = None,
                          all_annotation_images: List[IsotopeImages] = None):
    """
        Convert METASPACE datasets to AnnData object with optical image.

        Requires only a dataset (SMDataset), database (e.g. `("HMDB", "v4")`), and fdr.
        Output of `SMDataset.results` and `SMDataset.all_ion_images` can be provided to work,
        since the function does not allow to specify all parameters.

        This function return the AnnData object with optical image, and ion images transformed to fit optical image.

        Metadata is stored in `adata.uns`.
        The optical image is stored in `adata.uns["spatial"]["image"]["images"]["hires"]`

        :param ds: METASPACE dataset object
        :param database: Database to download results. Is not used if `results` and `all_annotation_images` is provided.
        :param fdr: FDR threshold. Is not used if `results` and `all_annotation_images` is provided.
        :param results: Output of `SMDataset.results()` call.
        :param all_annotation_images: Output of `SMDataset.all_ion_images()` call.
        :return: AnnData object including images, metadata, and annotation information.
        :raises ValueError: If no optical image is available for dataset.
        """

    if ds._gqclient.getRawOpticalImage(ds.id)['rawOpticalImage'] is None:
        raise ValueError('Dataset does not have an optical image')

    if results is None:
        res = ds.results(database=database, fdr=fdr)
    else:
        res = results.copy()

    res = m_results_replace_index(res)

    if all_annotation_images is None:
        ann_images = ds.all_annotation_images(fdr=fdr,
                                              database=database,
                                              only_first_isotope=True)
    else:
        ann_images = all_annotation_images.copy()

    # Sort them according to res
    img_dict = {image.formula + image.adduct: image for image in ann_images}

    im_array = np.array([img_dict[ion] for ion in res.index])

    # Remove additional dimension with placeholder for isotopes
    im_array = reshape_im_array(im_array)

    max_array = im_array.max(axis=1).max(axis=1)

    # Get optical images
    oi = ds.optical_images()

    # Replace im_array with transformed images
    im_array = np.array([np.asarray(oi.ion_image_to_optical(im)).sum(axis=2) for im in im_array])

    # Transform to X,Y representation
    xdata = im_array.reshape((im_array.shape[0], -1)).transpose()

    xy_coord, xy_coord_df = get_xy_coordinates(im_array, xdata)

    # Scaling data back to original
    xdata = (xdata / 3 / 255) * max_array

    adata = AnnData(X=xdata,
                    var=res,
                    obs=xy_coord_df,
                    obsm={'spatial': xy_coord})

    # Setting the BG image
    spatial_key = "spatial"
    library_id = "image"
    adata.uns[spatial_key] = {library_id: {}}
    adata.uns[spatial_key][library_id]["images"] = {}
    adata.uns[spatial_key][library_id]["images"] = {"hires": oi[0]}
    adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 1, "spot_diameter_fullres": 0.5}

    # Get METASPACE infos
    add_metadata(adata, ds)

    return adata
