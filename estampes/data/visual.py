"""Store data/constants for visualization.

The module stores data and/or constants relevant to different submodules
of the `visual` module.

Attributes
----------
BONDDATA
    Bond data for visualization.
MOLCOLS
    List of molecular colors with a sufficiently good contrast.
PATH_OBJ3D
    Path to files describing 3D objects and materials.
RAD_VIS_SCL
    Default scaling factors of radius for visualization.
"""

import os

# ================
# Module Constants
# ================

BONDDATA = {
    'rvis': 0.15,
    'rgb': (200, 200, 200)
}

MOLCOLS = [
    (0x1F, 0x77, 0xB4),
    (0xFF, 0x7F, 0x0C),
    (0x2C, 0xA0, 0x2C),
    (0xD6, 0x27, 0x28),
    (0x94, 0x67, 0xBD),
    (0x8C, 0x56, 0x4B),
    (0xE3, 0x77, 0xC2),
    (0x7F, 0x7F, 0x7F),
    (0xBC, 0xBD, 0x22),
    (0x17, 0xBE, 0xCF),
]

PATH_OBJ3D = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '3d_models')

RAD_VIS_SCL = .6
