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
MATERIALS
    Supported materials for visualization and software-specific names.
MODELS
    Representation models.
"""

import os

from estampes.data.colors import MPL_COLORS

# ================
# Module Constants
# ================

BONDDATA = {
    'rvis': 0.15,
    'rgb': (200, 200, 200)
}

MOLCOLS = list(MPL_COLORS.values())


PATH_OBJ3D = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '3d_models')

RAD_VIS_SCL = .6

MATERIALS = {
    'glass': {
        'povray': 'Mater_Glass'
    },
    'metal': {
        'povray': 'Mater_Metal'
    },
    'plastic': {
        'povray': 'Mater_Metal'
    },
}

MODELS = {
    'mol': {
        'balls': {
            'alias': ('ball', 'balls', 'ballsandsticks')
        },
        'spheres': {
            'alias': ('sphere', 'spheres', 'vdW')
        },
        'sticks': {
            'alias': ('stick', 'sticks')
        },
    },
    'vib': {
        'arrows': {
            'alias': ('arrow', 'arrows')
        },
        'spheres': {
            'alias': ('sphere', 'spheres')
        },
        'midarrows': {
            'alias': ('carrows', 'carrow', 'midarrows', 'midarrow')
        },
        'dualarrows': {
            'alias': ('2arrows', '2arrow', 'dualarrows', 'dualarrow')
        }
    }
}
