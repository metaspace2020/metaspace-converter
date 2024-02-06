from dataclasses import dataclass
from typing import Final

SPATIAL_KEY: Final = "spatial"
METASPACE_KEY: Final = "metaspace"
OBS_INDEX_NAME: Final = "ion_image_pixel"
VAR_INDEX_NAME: Final = "formula_adduct"
Shape2d = tuple[int, int]


@dataclass
class COL:
    ion_image_shape_y: Final = "ion_image_shape_y"
    ion_image_shape_x: Final = "ion_image_shape_x"
    ion_image_pixel_y: Final = "ion_image_pixel_y"
    ion_image_pixel_x: Final = "ion_image_pixel_x"


COORD_SYS_GLOBAL: Final = "global"
MICROMETER: Final = "micrometer"
OPTICAL_IMAGE_KEY: Final = "optical_image"
POINTS_KEY: Final = "maldi_points"
REGION_KEY: Final = "region"
INSTANCE_KEY: Final = "instance_id"
X: Final = "x"
Y: Final = "y"
C: Final = "c"
XY: Final = (X, Y)
YX: Final = (Y, X)
YXC: Final = (Y, X, C)

COLOCALIZATION = "colocalization"
