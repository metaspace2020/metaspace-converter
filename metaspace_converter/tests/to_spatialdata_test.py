import numpy as np
import pytest
from multiscale_spatial_image import MultiscaleSpatialImage

from metaspace_converter.constants import COORD_SYS_GLOBAL, OPTICAL_IMAGE_KEY, POINTS_KEY, X, Y

# The import metaspace_credentials is need by the fixture sm.
from metaspace_converter.tests.to_anndata_test import metaspace_credentials, sm  # noqa
from metaspace_converter.to_anndata import COL, METASPACE_KEY, get_ion_image_shape
from metaspace_converter.to_spatialdata import metaspace_to_spatialdata


@pytest.mark.parametrize(
    (
        "dataset_id",
        "database",
        "fdr",
        "add_optical_image",
        "optical_name_added",
        "add_points",
        "points_name_added",
    ),
    [
        ("2021-09-03_11h43m13s", ("CoreMetabolome", "v3"), 0.1, True, "optical", True, "maldi"),
        ("2021-09-03_11h43m13s", ("CoreMetabolome", "v3"), 0.1, False, None, False, None),
    ],
)
def test_metaspace_to_spatialdata(
    dataset_id,
    database,
    fdr,
    sm,
    add_optical_image,
    optical_name_added,
    add_points,
    points_name_added,
    tmp_path,
):
    actual = metaspace_to_spatialdata(
        dataset_id=dataset_id,
        database=database,
        fdr=fdr,
        sm=sm,
        add_optical_image=add_optical_image,
        optical_name_added=optical_name_added,
        add_points=add_points,
        points_name_added=points_name_added,
    )

    # Check table
    assert actual.table is not None
    dataset = sm.dataset(id=dataset_id)
    assert actual.table.n_obs == np.prod(get_ion_image_shape(dataset))
    assert actual.table.n_vars == len(dataset.annotations(fdr=fdr, database=database))
    assert {
        COL.ion_image_shape_y,
        COL.ion_image_shape_x,
        COL.ion_image_pixel_y,
        COL.ion_image_pixel_x,
    }.issubset(actual.table.obs.columns)
    assert np.all(actual.table.var["fdr"] <= fdr)
    actual_database = actual.table.uns[METASPACE_KEY]["metaspace_to_anndata_parameters"]["database"]
    assert actual_database == list(database)

    # Check optical image, if requested
    assert (optical_name_added in actual.images) == add_optical_image

    # Check points, if requested
    assert (points_name_added in actual.points) == add_points

    # Check transformations
    if add_optical_image and add_points:
        # All points should be positioned within the image bounds.
        points_over_image = actual.transform_element_to_coordinate_system(
            actual.points[points_name_added], target_coordinate_system=COORD_SYS_GLOBAL
        )[[Y, X]].compute()
        image_min = (0, 0)
        image_max = _get_image_shape(actual.images[optical_name_added])
        assert np.all(image_min <= points_over_image) and np.all(points_over_image <= image_max)

    # Check that it can be successfully written.
    actual.write(tmp_path / "sdata.zarr")


def _get_image_shape(image: MultiscaleSpatialImage) -> tuple[int, int]:
    assert isinstance(image, MultiscaleSpatialImage)
    sizes = image["scale0"].sizes
    return sizes["y"], sizes["x"]


@pytest.mark.parametrize(
    ("dataset_id", "database", "fdr"), [("2021-09-03_11h43m13s", ("CoreMetabolome", "v3"), 0.1)]
)
def test_metaspace_with_napari_spatialdata(dataset_id, database, fdr, sm, request):
    # A demo of viewing SpatialData with Napari.
    if (
        len(request.session.items) != 1
        or request.session.items[0].originalname != "test_metaspace_with_napari_spatialdata"
    ):
        pytest.skip("Manual test case. To run it, select it individually.")

    actual = metaspace_to_spatialdata(dataset_id=dataset_id, database=database, fdr=fdr)

    from napari_spatialdata import Interactive, QtAdataViewWidget

    interactive = Interactive(sdata=actual, headless=True)
    sdata_widget = interactive._sdata_widget
    anndata_widget = [
        w
        for w in interactive._viewer.window._dock_widgets["View (napari-spatialdata)"].children()
        if isinstance(w, QtAdataViewWidget)
    ][0]
    # 1. Select coordinate system
    sdata_widget.coordinate_system_widget._select_coord_sys(COORD_SYS_GLOBAL)
    sdata_widget.elements_widget._onItemChange(COORD_SYS_GLOBAL)
    # 2. Select image element
    sdata_widget._onClick(OPTICAL_IMAGE_KEY)
    # 3. Select points element
    sdata_widget._onClick(POINTS_KEY)
    # 4. Select observation for points
    anndata_widget.obs_widget._onAction([actual.table.var.index[0]])
    # Run Napari (blocking)
    interactive.run()
