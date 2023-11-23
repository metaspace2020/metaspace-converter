import json
import os
from pathlib import Path

import numpy as np
import pytest
from metaspace import SMInstance

from metaspace_converter.constants import COL, METASPACE_KEY, SPATIAL_KEY
from metaspace_converter.to_anndata import get_ion_image_shape, metaspace_to_anndata

METASPACE_DEFAULT_CONFIG_FILE = "~/.metaspace"
METASPACE_EMAIL_ENV = "METASPACE_EMAIL"
METASPACE_API_KEY_ENV = "METASPACE_API_KEY"


@pytest.fixture(scope="session")
def metaspace_credentials() -> dict[str, str]:
    if Path(METASPACE_DEFAULT_CONFIG_FILE).expanduser().exists():
        # Credentials are available. No need to read them, SMInstance will read the file.
        return {}
    else:
        return {
            "email": os.environ.get(METASPACE_EMAIL_ENV),
            "api_key": os.environ.get(METASPACE_API_KEY_ENV),
        }


@pytest.fixture(scope="session")
def sm(metaspace_credentials) -> SMInstance:
    return SMInstance(
        email=metaspace_credentials.get("email"), api_key=metaspace_credentials.get("api_key")
    )


@pytest.mark.parametrize(
    ("dataset_id", "database", "fdr", "metadata_as_obs", "add_optical_image"),
    [
        ("2021-09-03_11h43m13s", ("CoreMetabolome", "v3"), 0.1, False, False),
        ("2021-09-03_11h43m13s", ("CoreMetabolome", "v3"), 0.1, True, False),
        ("2021-09-03_11h43m13s", ("CoreMetabolome", "v3"), 0.1, False, True),
    ],
)
def test_metaspace_to_anndata(
    dataset_id, database, fdr, sm, metadata_as_obs, add_optical_image, tmp_path
):
    actual = metaspace_to_anndata(
        dataset_id=dataset_id,
        database=database,
        fdr=fdr,
        sm=sm,
        metadata_as_obs=metadata_as_obs,
        add_optical_image=add_optical_image,
    )

    dataset = sm.dataset(id=dataset_id)
    assert actual.n_obs == np.prod(get_ion_image_shape(dataset))
    assert actual.n_vars == len(dataset.annotations(fdr=fdr, database=database))
    assert {
        COL.ion_image_shape_y,
        COL.ion_image_shape_x,
        COL.ion_image_pixel_y,
        COL.ion_image_pixel_x,
    }.issubset(actual.obs.columns)
    assert np.all(actual.var["fdr"] <= fdr)
    if metadata_as_obs:
        assert all(
            [
                json.loads(row["metaspace.metaspace_to_anndata_parameters"])["database"]
                == list(database)
                for _, row in actual.obs.iterrows()
            ]
        )
    else:
        actual_database = actual.uns[METASPACE_KEY]["metaspace_to_anndata_parameters"]["database"]
        assert actual_database == list(database)
    if add_optical_image:
        assert SPATIAL_KEY in actual.obsm
        assert SPATIAL_KEY in actual.uns
        assert isinstance(actual.uns[SPATIAL_KEY]["image"]["images"]["hires"], np.ndarray)

    # Check that it can be successfully written.
    actual.write_h5ad(tmp_path / "adata.h5ad")
    actual.write_zarr(tmp_path / "adata.zarr")
