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

from pyabsorp.air import air_properties

from pyabsorp.models import delany_bazley, delany_bazley_absorption, \
    rayleigh, rayleigh_absorption, biot_allard, biot_allard_absorption, \
    johnson_champoux, johnson_champoux_absorption

__author__ = "Michael Markus Ackermann"
__date__ = "27 March 2020"
__version__ = "1.0.0"

# package submodules and scripts to be called as pyabsorp.something
__all__ = [  # Functions
            'air_properties',
            'delany_bazley',
            'delany_bazley_absorption',
            'rayleigh',
            'rayleigh_absorption',
            'biot_allard',
            'biot_allard_absorption',
            'johnson_champoux',
            'johnson_champoux_absorption']
