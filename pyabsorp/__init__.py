"""
PYABSORP
=====================================================
Author:
    Michael Markus Ackermann - dev.toktom@outlook.com
PyAbsorp:
    This is a package developed to be use to find the Sound Absorption
    Coefficient through some implemented models, like  Biot-Allard,
    Johnson-Champoux and others.
    In order to provide such functionalities we require a few packages
    to be installed:
    - Numpy
    - Scipy

    You can find out everything available reading the submodules documentation


For further information, check the specific module, class, method or function
documentation.
"""

from pyabsorp.air import air_properties, AirProperties
from pyabsorp.absorption import absorption_coefficient
from pyabsorp.version import __author__, __date__, __version__
from pyabsorp.models import delany_bazley, rayleigh, biot_allard, johnson_champoux
from pyabsorp.material import Material

# Just to prevent "unused" warnings.
assert __author__ and __date__ and __version__

# package submodules and scripts to be called as pyabsorp.something
__all__ = [  # Functions
            'absorption_coefficient',
            'air_properties',
            'delany_bazley',
            'rayleigh',
            'biot_allard',
            'johnson_champoux',

            # Classes
            'AirProperties',
            'Material']
