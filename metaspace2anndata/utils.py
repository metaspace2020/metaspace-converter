import numpy as np
import pandas as pd
from metaspace.sm_annotation_utils import SMDataset
from anndata import AnnData

def reshape_im_array(im_array: np.ndarray):
    return im_array.reshape((im_array.shape[0],
                             im_array.shape[2],
                             im_array.shape[3]))


def get_xy_coordinates(im_array: np.ndarray,
                       xdata: np.ndarray):
    xy_coord = np.indices(im_array.shape)[1:, :].reshape((2, -1))
    xy_coord = xy_coord[:, range(xdata.shape[0])].transpose()
    xy_coord = xy_coord[:, [1, 0]]

    xy_coord_df = pd.DataFrame(xy_coord, columns=['x', 'y'])

    return xy_coord, xy_coord_df


def m_results_replace_index(result: pd.DataFrame):
    res = result.reset_index()
    res['formula_adduct'] = res.reset_index()['formula'].str.cat(res.reset_index()['adduct'])
    return res.set_index('formula_adduct')


def add_metadata(adata: AnnData, ds: SMDataset):

    adata.uns['metaspace_metadata'] = ds.metadata
    adata.uns['metaspace_id'] = ds.id
    adata.uns['metaspace_database_details'] = ds.database_details
    adata.uns['metaspace_config'] = ds.config
    adata.uns['metaspace_name'] = ds.name
    adata.uns['metaspace_adducts'] = ds.adducts
    adata.uns['metaspace_databases'] = ds.databases
    adata.uns['metaspace_group'] = ds.group
    adata.uns['metaspace_polarity'] = ds.polarity
    adata.uns['metaspace_principal_investigator'] = ds.principal_investigator
    adata.uns['metaspace_projects'] = ds.projects
    adata.uns['metaspace_status'] = ds.status
    adata.uns['metaspace_submitter'] = ds.submitter


