import dataclasses
import json
from collections.abc import Collection, Iterator, Mapping, Sequence
from dataclasses import is_dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd
from metaspace.sm_annotation_utils import SMDataset


def empty_columns_to_nan(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts all empty columns from None (dtype object) to NaN (dtype float).

    Args:
        df: A dataframe

    Returns:
        The modified dataframe. The original dataframe is modified, not a copy.
    """
    for c in df.columns:
        if pd.api.types.infer_dtype(df[c]) == "empty":
            df[c] = np.nan
    return df


def has_optical_image(dataset: SMDataset) -> bool:
    """
    Check whether a METASPACE dataset has an optical image
    """
    try:
        optical_images = dataset.optical_images()
    except TypeError:
        return False
    return optical_images is not None and len(optical_images) > 0


def is_non_str_sequence(obj: Any) -> bool:
    """
    Returns whether an object is a sequence (list and similar), excluding strings.

    Args:
        obj: The object to check.

    Returns:
        True if the object is a sequence.
    """
    return isinstance(obj, Sequence) and not isinstance(obj, str)


def iter_property_paths(
    obj: Any,
    delimiter: str = "/",
    include_none: bool = True,
    dont_traverse_types: tuple[type, ...] = tuple(),
    dont_descend_paths: Collection[str] = tuple(),
    _path: str = "",
) -> Iterator[tuple[str, Any]]:
    """
    Iterates recursively over properties of nested objects.

    Args:
        obj: A class instance, mapping or sequence.
        delimiter: The delimiter to use. Should not be contained in any name.
        include_none: Whether to yield attributes which have None assigned, or skip them.
        dont_traverse_types: Iterable data types which should not be traversed any deeper.
            These are yielded as values.
        dont_descend_paths: Paths with iterable objects that should be preserved.

    Returns:
        An iterator for iterating over paths and objects found at each path. It can also be
        converted to a dictionary.
    """
    prefix = _path if _path == "" else f"{_path}{delimiter}"
    if dont_traverse_types and isinstance(obj, dont_traverse_types):
        if obj is None and not include_none:
            return
        yield _path, obj
    elif _path in dont_descend_paths:
        if obj is None and not include_none:
            return
        yield _path, obj
    # Traversable types: list, dict, dataclasses/pydantic
    elif is_non_str_sequence(obj):
        for i, value in enumerate(obj):
            yield from iter_property_paths(
                value,
                delimiter=delimiter,
                include_none=include_none,
                dont_traverse_types=dont_traverse_types,
                dont_descend_paths=dont_descend_paths,
                _path=f"{prefix}{i}",
            )
    elif isinstance(obj, Mapping):
        for key, value in obj.items():
            yield from iter_property_paths(
                value,
                delimiter=delimiter,
                include_none=include_none,
                dont_traverse_types=dont_traverse_types,
                dont_descend_paths=dont_descend_paths,
                _path=f"{prefix}{key}",
            )
    elif is_dataclass(obj):
        keys = (
            [field.name for field in dataclasses.fields(obj)]
            if is_dataclass(obj)
            else obj.__fields__.keys()
        )
        for key in keys:
            # Object with attributes
            value = getattr(obj, key)
            yield from iter_property_paths(
                value,
                delimiter=delimiter,
                include_none=include_none,
                dont_traverse_types=dont_traverse_types,
                dont_descend_paths=dont_descend_paths,
                _path=f"{prefix}{key}",
            )
    # Non-traversable/primitive types
    else:
        if obj is None and not include_none:
            # Don't yield values that are not set.
            return
        yield _path, obj


def stringify_list_columns(
    df: pd.DataFrame, columns: Iterable[str] = tuple(), all: bool = False
) -> pd.DataFrame:
    """
    If the dataframe column is a list, stringify it into JSON.

    Args:
        df: A dataframe with a column that contains lists (e.g. of molecule names)
        columns: The names of the columns containing lists
        all: Whether to check and stringify all columns containing list-like objects

    Returns:
        A dataframe where all values in specified columns are JSON stringified.
    """
    if not df.empty:
        df = df.copy()
        if all:
            columns = df.columns
        for column in columns:
            if column in df.columns and pd.api.types.is_list_like(df[column].iloc[0]):
                df[column] = df[column].apply(lambda l: json.dumps(l))
    return df


def transform_coords(coords: np.ndarray, transform_matrix_2d: np.ndarray) -> np.ndarray:
    """
    Transform a list of 2D coordinates using a 3×3 transformation matrix.

    Like in scikit-image, but without adding that dependency.

    See skimage.transform._geometric.ProjectiveTransform.__call__

    Args:
        coords: An array of coordinates with shape (number of points, 2)
        transform_matrix_2d: A transformation matrix of shape (3, 3)

    Returns:
        The list of transformed coordinates as Numpy array
    """
    # From skimage.transform._geometric
    # Convert to homogeneous coordinates
    src_h = np.hstack([coords, np.ones((len(coords), 1))])
    dst_h = src_h @ transform_matrix_2d.T
    # Below, we will divide by the last dimension of the homogeneous
    # coordinate matrix. In order to avoid division by zero,
    # we replace exact zeros in this column with a very small number.
    dst_h[dst_h[:, -1] == 0, -1] = np.finfo(float).eps
    # rescale to homogeneous coordinates
    dst = dst_h[:, :-1] / dst_h[:, -1:]
    return dst
