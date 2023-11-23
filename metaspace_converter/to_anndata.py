from typing import Container, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from metaspace import SMInstance
from metaspace.sm_annotation_utils import (
    DEFAULT_DATABASE,
    DatasetNotFound,
    IsotopeImages,
    SMDataset,
)

from metaspace_converter.constants import (
    COL,
    METASPACE_KEY,
    OBS_INDEX_NAME,
    SPATIAL_KEY,
    VAR_INDEX_NAME,
    Shape2d,
)
from metaspace_converter.utils import (
    empty_columns_to_nan,
    has_optical_image,
    iter_property_paths,
    stringify_list_columns,
    transform_coords,
)


def metaspace_to_anndata(
    dataset: Optional[SMDataset] = None,
    dataset_id: Optional[str] = None,
    database: tuple[str, str] = None,
    fdr: Optional[float] = 0.1,
    use_tic: bool = False,
    metadata_as_obs: bool = False,
    add_optical_image: bool = False,
    sm: Optional[SMInstance] = None,
    **annotation_filter,
) -> AnnData:
    """
    Downloads a METASPACE dataset to an :py:class:`AnnData` object.

    See: https://metaspace2020.eu/about

    Args:
        dataset: A METASPACE dataset instance. If not provided, the `dataset_id` must be given.
        dataset_id: The unique ID of a dataset on METASPACE, e.g. "2021-09-03_11h43m13s"
        database: A single METASPACE database given as a tuple of name and version. Usually it is
            displayed on METASPACE as "HMDB – v4" which corresponds to `("HMDB", "v4")`.
        fdr: Returns only annotations for which the false discovery rate is less or equal to this
            limit.
        use_tic: When True, the output values will be scaled by the total ion count per pixel and
            will be in 0.0 to 1.0 range.
        metadata_as_obs: Whether to store metadata in the `obs` dataframe instead of `uns`. For a
            single METASPACE dataset, metadata is the same for all pixels, so it would be duplicated
            for all `obs`. When combining multiple datasets, it would be preserved in `obs` but
            not in `uns`.
        add_optical_image: Whether to embed the optical image for SquidPy
        sm: Optionally a cached SMInstance
        annotation_filter: Additional keyword arguments passed to the METASPACE API.

    Returns:
        An :py:class:`AnnData` object with
            * ion intensities: ``.X``
            * ion image pixel coordinates: ``.obs[["ion_image_pixel_x", "ion_image_pixel_y"]]``
            * spatial coordinates: ``.obsm["spatial"]``
            * ion properties: ``.var``, for example "formula", "adduct", "mz", "fdr",
              "moleculeNames", "moleculeIds", "intensity"…
            * METASPACE metadata: ``.uns["metaspace"]`` if not ``metadata_as_obs``
            * SquidPy metadata: ``.uns["spatial"]`` if ``add_optical_image``

    Raises:
        ValueError: If something is wrong with the input data or parameters
        metaspace.sm_annotation_utils.DatasetNotFound: If the dataset ID is not found
        ConnectionError: For any other network connectivity problems
    """
    # This must stay at the top above `locals()`, so that `parameters` includes no later defined
    # variables, but `dataset_id` definitely has a value.
    if dataset is not None:
        dataset_id = dataset.id
    # Capture the function parameters (except non-serializable ones) so we can include them
    # in the resulting AnnData.
    parameters = _dict_without_keys(locals(), keys=("sm", "dataset"))

    dataset = get_dataset(dataset=dataset, dataset_id=dataset_id, sm=sm)
    if dataset.status != "FINISHED":
        raise DatasetNotFound(f"Dataset with ID {dataset.id} has not finished processing.")

    if database is None:
        database = DEFAULT_DATABASE

    # Download annotations
    annotations = dataset.results(database=database, fdr=fdr, **annotation_filter)
    annotations = _add_annotations_index(annotations, index_name=VAR_INDEX_NAME)
    annotations = _normalize_annotations_for_serialization(annotations)

    # Download ion images
    isotope_images = dataset.all_annotation_images(
        only_first_isotope=True,
        fdr=fdr,
        database=database,
        scale_intensity="TIC" if use_tic else True,
        **annotation_filter,
    )
    if not isotope_images:
        raise ValueError(
            f"No isotope images available for dataset {dataset.id} and database "
            f"{database[0]} – {database[1]}. Was the database selected for processing on METASPACE?"
        )
    assert len(annotations) == len(isotope_images)

    # Sort them matching the annotations.
    isotope_images = _sort_isotope_images_like(isotope_images, annotations.index)

    # Create X matrix (all ion pixels flattened to primary axis)
    shape = get_ion_image_shape(dataset)
    ion_intensities = _create_anndata_x(isotope_images, shape)

    # Create observations dataframe
    obs = _create_obs_with_coordinates(shape)

    # Scanpy looks up (x, y) coordinates from obsm.
    obsm = {SPATIAL_KEY: obs[[COL.ion_image_pixel_x, COL.ion_image_pixel_y]].values}

    adata = AnnData(X=ion_intensities, var=annotations, obs=obs, obsm=obsm)

    # Add metadata
    if not metadata_as_obs:
        adata.uns[METASPACE_KEY] = _create_uns_for_metaspace(dataset, parameters=parameters)
    else:
        metadata_df = _create_obs_for_metaspace(
            dataset, obs_names=adata.obs.index, parameters=parameters
        )
        adata.obs = pd.concat([obs, metadata_df], axis=1)

    # Add optical image
    if add_optical_image and has_optical_image(dataset):
        _add_obsm_and_uns_for_squidpy(
            adata,
            dataset=dataset,
            pixel_coords_xy=obs[[COL.ion_image_pixel_x, COL.ion_image_pixel_y]],
        )

    return adata


def create_annotation_id(formula: str, adduct: str) -> str:
    return f"{formula}{adduct}"


def _add_annotations_index(df: pd.DataFrame, index_name: str = VAR_INDEX_NAME) -> pd.DataFrame:
    df = df.reset_index()
    df[index_name] = df.apply(lambda row: create_annotation_id(row.formula, row.adduct), axis=1)
    return df.set_index(index_name)


def _normalize_annotations_for_serialization(df: pd.DataFrame) -> pd.DataFrame:
    # AnnData cannot save dataframes containing Python objects (like lists).
    # As workaround, we stringify to JSON.
    df = stringify_list_columns(df, columns=["isotopeImages", "moleculeNames", "moleculeIds"])
    # AnnData can also not save columns containing only None (dtype object, tries to save as str)
    # "Can't implicitly convert non-string objects to strings"
    return empty_columns_to_nan(df)


def get_dataset(
    dataset: Optional[SMDataset] = None,
    dataset_id: Optional[str] = None,
    sm: Optional[SMInstance] = None,
) -> SMDataset:
    """
    Fetch a dataset instance, if not provided
    """
    if dataset is not None:
        return dataset
    else:
        if dataset_id is None:
            raise ValueError("Either `dataset` or `dataset_id` must be provided")
        if sm is None:
            sm = SMInstance()
        return sm.dataset(id=dataset_id)


def get_ion_image_shape(
    dataset: Optional[SMDataset] = None,
    dataset_id: Optional[str] = None,
    sm: Optional[SMInstance] = None,
) -> Shape2d:
    """
    Get the shape of a METASPACE dataset (in Numpy order).

    Args:
        dataset: A METASPACE dataset instance. If not provided, a dataset ID must be provided.
        dataset_id: An optional METASPACE dataset ID, if no dataset instance is provided.
        sm: An optional SMInstance, otherwise one will be created.
    """
    dataset = get_dataset(dataset=dataset, dataset_id=dataset_id, sm=sm)
    return dataset.image_size["y"], dataset.image_size["x"]


def _sort_isotope_images_like(
    isotope_images: list[IsotopeImages], index: pd.Index
) -> list[IsotopeImages]:
    images_dict = {}
    for isotope_image in isotope_images:
        annotation_id = create_annotation_id(isotope_image.formula, isotope_image.adduct)
        images_dict[annotation_id] = isotope_image
    # Return them in the requested order.
    return [images_dict[key] for key in index]


def _create_anndata_x(isotope_images: list[IsotopeImages], shape: Shape2d) -> np.ndarray:
    n_pixels = np.prod(shape)
    n_ions = len(isotope_images)
    y_x_ion_stack = np.dstack([_isotope_image_to_array(img, shape=shape) for img in isotope_images])
    return y_x_ion_stack.reshape(n_pixels, n_ions)


def _isotope_image_to_array(ion_image: IsotopeImages, shape: Shape2d) -> np.ndarray:
    if len(ion_image) == 0 or ion_image[0] is None:
        return np.zeros(shape=shape)
    else:
        return ion_image[0]


def _create_obs_with_coordinates(shape: Shape2d) -> pd.DataFrame:
    n_pixels = np.prod(shape)
    pixel_coords_yx = all_image_pixel_coordinates(shape)
    return pd.DataFrame(
        {
            COL.ion_image_pixel_x: pixel_coords_yx[:, 1],
            COL.ion_image_pixel_y: pixel_coords_yx[:, 0],
            COL.ion_image_shape_y: shape[0],
            COL.ion_image_shape_x: shape[1],
        },
        index=pd.RangeIndex(n_pixels, name=OBS_INDEX_NAME).astype(str),
    )


def all_image_pixel_coordinates(shape: Shape2d) -> np.ndarray:
    """
    Returns a list of Y, X coordinates (Numpy order) for every pixel of an image.

    Args:
        shape: The shape of a 2D image

    Returns:
        Array of shape (number of pixels, 2)
    """
    ys_and_xs = np.indices(shape)
    return np.dstack(ys_and_xs).reshape((-1, len(shape)))


def _create_default_spatial_uns(
    image: Optional[np.ndarray] = None,
    scale_factor: float = 1.0,
    spot_diameter: float = 1,
    library_id: str = "image",
    image_key: str = "hires",
) -> dict:
    return {
        library_id: {
            "images": {image_key: image} if image is not None else {},
            "scalefactors": {
                f"tissue_{image_key}_scalef": scale_factor,
                "spot_diameter_fullres": spot_diameter,
            },
        }
    }


def _create_uns_for_metaspace(dataset: SMDataset, parameters: Optional[dict] = None) -> dict:
    # Database details contain instances of MolecularDB, which are not JSON-serializable.
    # Also, AnnData does  not support lists of dicts, so we have to convert to dict of dicts.
    # See issue https://github.com/scverse/anndata/issues/708
    plain_db_details = {str(idx): db._info for idx, db in enumerate(dataset.database_details)}
    plain_projects = {str(idx): project for idx, project in enumerate(dataset.projects)}
    plain_parameters = {k: list(v) if isinstance(v, tuple) else v for k, v in parameters.items()}
    return {
        "id": dataset.id,
        "name": dataset.name,
        "adducts": dataset.adducts,
        "databases": dataset.databases,
        "database_details": plain_db_details,
        "polarity": dataset.polarity,
        "image_size": dataset.image_size,
        "submitter": dataset.submitter,
        "group": dataset.group,
        "principal_investigator": dataset.principal_investigator,
        "projects": plain_projects,
        "status": dataset.status,
        "metadata": dict(dataset.metadata),
        "config": dataset.config,
        # Include parameters used to generate this AnnData
        f"{metaspace_to_anndata.__name__}_parameters": plain_parameters,
    }


def _create_obs_for_metaspace(
    dataset: SMDataset, obs_names: pd.Index, parameters: Optional[dict] = None
) -> pd.DataFrame:
    metadata = _create_uns_for_metaspace(dataset=dataset, parameters=parameters)
    flattened_metadata = iter_property_paths(
        obj={METASPACE_KEY: metadata},
        delimiter=".",
        dont_traverse_types=(tuple, list),
        dont_descend_paths=[f"{METASPACE_KEY}.{metaspace_to_anndata.__name__}_parameters"],
    )
    metadata_series = pd.Series(dict(flattened_metadata))
    obs = pd.DataFrame([metadata_series] * len(obs_names), index=obs_names)
    # METASPACE metadata contains a lot of list values.
    # AnnData cannot save dataframes containing Python objects (like lists).
    # As workaround, we stringify to JSON.
    obs = stringify_list_columns(obs, all=True)
    return empty_columns_to_nan(obs)


def _add_obsm_and_uns_for_squidpy(adata: AnnData, dataset: SMDataset, pixel_coords_xy: np.ndarray):
    optical_images = dataset.optical_images()
    image_yx_rgb = optical_images[0]
    # Transformation from ion image coordinate system to optical image (x, y order!).
    matrix_xy = optical_images._transforms[0]
    coords_on_optical_xy = transform_coords(pixel_coords_xy, matrix_xy)
    # SquidPy looks up (x, y) coordinates from obsm.
    adata.obsm[SPATIAL_KEY] = coords_on_optical_xy
    # SquidPy uses spatial information from uns.
    adata.uns[SPATIAL_KEY] = _create_default_spatial_uns(image=image_yx_rgb, spot_diameter=0.5)


def _dict_without_keys(dictionary: dict, keys: Container) -> dict:
    # Exclude keys from a dictionary
    return {k: v for k, v in dictionary.items() if k not in keys}
