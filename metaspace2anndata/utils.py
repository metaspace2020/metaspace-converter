import numpy as np
import pandas as pd


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
