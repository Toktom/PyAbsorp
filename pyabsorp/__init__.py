# -*- coding: utf-8 -*-
"""
PYABSORP
========

Author:
    Michael Markus Ackermann - dev.toktom@outlook.com
PyAbsorp:
    This is a package developed to be use to find the Sound Absorption
    Coefficient through some implemented models, like Biot-Allard,
    Johnson-Champoux and others.
    In order to provide such functionalities we require a few packages
    to be installed:
    - Numpy
    - Scipy

    You can find out everything available reading the submodules documentation


For further information, check the specific module, class, method or function
documentation.
"""

from pyabsorp.version import __author__, __date__, __version__
from pyabsorp.models import delany_bazley, rayleigh, biot_allard, johnson_champoux
from pyabsorp.classes import Air, Material
from pyabsorp.functions import air_properties, air_impedance, air_density, \
    specific_heat_constant_pressure, specific_heat_constant_volume,\
    specific_heat_ratio, sound_speed, ceilsius_to_kelvin, pierce, \
    prandtl, viscosity, absorption_coefficient
from pyabsorp.utils import load_object_with_pickle, save_object_with_pickle

# Just to prevent "unused" warnings.
assert __author__ and __date__ and __version__

__all__ = [
    # Functions
    'absorption_coefficient',
    'delany_bazley',
    'rayleigh',
    'biot_allard',
    'johnson_champoux',
    'air_properties',
    'air_impedance',
    'air_density',
    'specific_heat_constant_pressure',
    'specific_heat_constant_volume',
    'specific_heat_ratio',
    'sound_speed',
    'ceilsius_to_kelvin',
    'pierce',
    'prandtl',
    'viscosity',
    'load_object_with_pickle',
    'save_object_with_pickle',
    # Classes
    'Material',
    'Air']
