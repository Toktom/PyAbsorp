"""
Author: Michael Markus Ackermann
================================

Models:
-------

    This submodule has some of the sound absorption coefficient models
    implemented, like the models from Delany-Bazley, Rayleigh, Biot-Allard
    and Johnson-Champoux, where Delany-Bazley and Johnson-Champoux have also
    the possibility to use their modified versions.

    Available functions:
    ---------------------

        >>> pyabsorp.delany_bazley(flow_resis, air_dens, ..., var='default')
        >>> pyabsorp.delany_bazley_absorption(zc, kc, d, z_air)
        >>> pyabsorp.rayleigh(flow_resis, air_dens, ..., freq)
        >>> pyabsorp.rayleigh_absorption(zc, kc, d, z_air)
        >>> pyabsorp.biot_allard(flow_resis, air_dens, ..., freq)
        >>> pyabsorp.biot_allard_absorption(zc, kc, d, z_air)
        >>> pyabsorp.johnson_champoux(flow_resis, air_dens, ... var='default')
        >>> pyabsorp.johnson_champoux_absorption(zc, kc, d, z_air)


    For further information, check the function specific documentation.
"""

from .delany_bazley import delany_bazley, delany_bazley_absorption
from .rayleigh import rayleigh, rayleigh_absorption
from .biot_allard import biot_allard, biot_allard_absorption
from .johnson_champoux import johnson_champoux, johnson_champoux_absorption

__all__ = [  # Functions
            'delany_bazley',
            'delany_bazley_absorption',
            'rayleigh',
            'rayleigh_absorption',
            'biot_allard',
            'biot_allard_absorption',
            'johnson_champoux',
            'johnson_champoux_absorption']
