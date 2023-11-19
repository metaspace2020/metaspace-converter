from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import pytest

from metaspace_converter.utils import (
    empty_columns_to_nan,
    iter_property_paths,
    stringify_list_columns,
    transform_coords,
)


def test_empty_columns_to_nan():
    df = pd.DataFrame({c: [v] * 2 for c, v in zip("abcdef", [1.0, 0.0, 0, object(), "", None])})
    df.loc[1, "a"] = np.nan
    df_orig = df.copy()
    empty_columns_to_nan(df)
    assert df["f"].dtype == float
    assert np.all(pd.isna(df["f"]))
    unchanged_columns = df_orig.columns != "f"
    pd.testing.assert_frame_equal(df.loc[:, unchanged_columns], df_orig.loc[:, unchanged_columns])


@dataclass
class A:
    a: Optional[int] = None
    b: Optional[int] = None
    c: Optional[int] = None


@dataclass
class C:
    l: Optional[list[A]] = None
    d: Optional[dict[str, A]] = None


@pytest.mark.parametrize(
    ("obj", "include_none", "expected"),
    [
        ("a", True, [("", "a")]),
        ({"a": 1, "b": 2}, True, [("a", 1), ("b", 2)]),
        ([3, 2], True, [("0", 3), ("1", 2)]),
        (A(a=1, b=2), False, [("a", 1), ("b", 2)]),
        (A(a=1, b=2), True, [("a", 1), ("b", 2), ("c", None)]),
        (
            C(l=[A(a=1)], d={"a": 1, "b": A(a=1, b=0)}),
            False,
            [("l/0/a", 1), ("d/a", 1), ("d/b/a", 1), ("d/b/b", 0)],
        ),
    ],
)
def test_iter_property_paths(obj, include_none, expected):
    actual = list(iter_property_paths(obj, include_none=include_none))
    assert actual == expected


@pytest.mark.parametrize(
    ("obj", "delimiter", "expected"),
    [
        ("a", ".", [("", "a")]),
        (
            C(l=[A(a=1)], d={"a": 1, "b": A(a=1, b=0)}),
            ".",
            [("l.0.a", 1), ("d.a", 1), ("d.b.a", 1), ("d.b.b", 0)],
        ),
    ],
)
def test_iter_property_paths_delimiter(obj, delimiter, expected):
    actual = list(iter_property_paths(obj, delimiter=delimiter, include_none=False))
    assert actual == expected


@pytest.mark.parametrize(
    ("obj", "dont_traverse_types", "dont_descend_paths", "expected"),
    [
        (C(l=[A(a=1)], d={"b": A(a=1)}), (list,), [], [("l", [A(a=1)]), ("d/b/a", 1)]),
        (C(l=[A(a=1)], d={"b": A(a=1)}), (dict,), [], [("l/0/a", 1), ("d", {"b": A(a=1)})]),
        (C(l=[A(a=1)], d={"b": A(a=1)}), tuple(), ["d/b"], [("l/0/a", 1), ("d/b", A(a=1))]),
    ],
)
def test_iter_property_paths_exclusion(obj, dont_traverse_types, dont_descend_paths, expected):
    actual = list(
        iter_property_paths(
            obj,
            dont_traverse_types=dont_traverse_types,
            dont_descend_paths=dont_descend_paths,
            include_none=False,
        )
    )
    assert actual == expected


@pytest.mark.parametrize(
    ("df", "expected", "columns"),
    [
        (
            pd.DataFrame({"col1": [["db", "v1"]], "col2": [[{}]]}),
            pd.DataFrame({"col1": ['["db", "v1"]'], "col2": [[{}]]}),
            ["col1"],
        ),
        (
            pd.DataFrame({"col1": [["db", "v1"]], "col2": [[{}]]}),
            pd.DataFrame({"col1": ['["db", "v1"]'], "col2": ["[{}]"]}),
            ["col1", "col2"],
        ),
    ],
)
def test_stringify_list_columns(df, columns, expected):
    actual = stringify_list_columns(df, columns=columns)

    pd.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    ("coordinates", "matrix", "expected"),
    [
        ([[0, 0], [1, 2]], np.array([[1, 0, 4], [0, 1, 5], [0, 0, 1]]), [[4, 5], [5, 7]]),
        ([[0, 0], [1, 2]], np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]), [[0, 0], [-2, 1]]),
    ],
)
def test_transform_coords(coordinates, matrix, expected):
    actual = transform_coords(np.asarray(coordinates), matrix)

    np.testing.assert_allclose(actual, expected)
