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

        >>> pyabsorp.delany_bazley(flow_resis, air_dens, ..., var='default')#
        >>> pyabsorp.rayleigh(flow_resis, air_dens, ..., freq)
        >>> pyabsorp.biot_allard(flow_resis, air_dens, ..., freq)
        >>> pyabsorp.johnson_champoux(flow_resis, air_dens, ... var='default')


    For further information, check the function specific documentation.
"""

from .delany_bazley import delany_bazley
from .rayleigh import rayleigh
from .biot_allard import biot_allard
from .johnson_champoux import johnson_champoux

__all__ = [  # Functions
            'delany_bazley',
            'rayleigh',
            'biot_allard',
            'johnson_champoux']
