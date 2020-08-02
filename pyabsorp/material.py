#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide basic interface to handle a single material being studied.


Created on Wed Jul 29 23:09:54 2020

author: João Vitor Gutkoski Paes
"""

import numpy as np
from pyabsorp.air import AirProperties
from pyabsorp.absorption import absorption_coefficient
from pyabsorp.models import delany_bazley, rayleigh, biot_allard, johnson_champoux


class Material(object):
    """Basic material object interface."""

    def __init__(self, thick: float, freq: np.ndarray,
                 air: AirProperties, *, poros: float = None,
                 tortus: float = None, flowres: float = None,
                 thermlen: float = None, visclen: float = None,
                 shape: str = None, thermperm: float = None):
        """
        Representation of a material being studied.

        This class provides an interface to all models of absorption coefficient
        available in `PyAbsorp`.

        For now each parameter, except `air`, can be easily set by an assignment
        operation (a = b), but this should be improved to at least some basic type
        and numerical range checking.

        Most parameters are derived from laboratory tests and are described as LP,
        for Laboratory Parameter. This class provides interface to semi-empirical
        and analytical models, no regression based on measurements are available.
        At least not yet.

        All modelling is made on the frequency domain.


        Parameters
        ----------
        thick : float
            Material thickness.
        freq : np.ndarray
            Array of frequencies.
        air : AirProperties
            Air acoustical properties.
        poros : float, optional
            Open porosity, LP. The default is None.
        tortus : float, optional
            Tortuosity, LP. The default is None.
        flowres : float, optional
            Static flow resistivity, LP. The default is None.
        thermlen : float, optional
            Thermal characteristic length, LP. The default is None.
        visclen : float, optional
            Viscous characteristic length, LP. The default is None.
        shape : str, optional
            Shape of the pore, LP. The default is None.
        thermperm : float, optional
            Static thermal permeability, LP. The default is None.

        Returns
        -------
        None.

        """
        self._air = air
        self.thickness = thick
        self.frequencies = np.float32(freq)
        self.porosity = poros
        self.tortuosity = tortus
        self.flowResistivity = flowres
        self.thermalLength = thermlen
        self.viscousLength = visclen
        self.poreShape = shape
        self.thermalPerm = thermperm
        self._kc = self._zc = self._absorp = None
        return

    @property
    def thickness(self):
        return self._thick

    @thickness.setter
    def thickness(self, thick):
        self._thick = thick
        return

    @property
    def porosity(self):
        return self._poros

    @porosity.setter
    def porosity(self, poros):
        self._poros = poros
        return

    @property
    def poreShape(self):
        return self._shape

    @poreShape.setter
    def poreShape(self, shape):
        self._shape = shape
        return

    @property
    def tortuosity(self):
        return self._tortus

    @tortuosity.setter
    def tortuosity(self, tortus):
        self._tortus = tortus
        return

    @property
    def flowResistivity(self):
        return self._flowres

    @flowResistivity.setter
    def flowResistivity(self, flowres):
        self._flowres = flowres
        return

    @property
    def thermalLength(self):
        return self._thermlen

    @thermalLength.setter
    def thermalLength(self, thermlen):
        self._thermlen = thermlen
        return

    @property
    def thermalPerm(self):
        return self._thermperm

    @thermalPerm.setter
    def thermalPerm(self, thermperm):
        self._thermperm = thermperm
        return

    @property
    def viscousLength(self):
        return self._visclen

    @viscousLength.setter
    def viscousLength(self, visclen):
        self._visclen = visclen
        return

    @property
    def frequencies(self):
        return self._freq

    @frequencies.setter
    def frequencies(self, freq):
        self._freq = freq
        return

    @property
    def air(self):
        return self._air

    @property
    def impedance(self):
        return self._zc

    @property
    def waveNum(self):
        return self._kc

    @property
    def absorption(self):
        return self._absorp

    def estimate_absorption(self, method: str, var: str = 'default'):
        """
        Estimate material absorption based on `method`.

        The material will hold the resulting `absorption` coefficients and the
        respective characteristic `impedance` and wave number (`waveNum`).

        Only the result of one call can be held. This means that a comparison between
        methods, or method variations must save separatedly each `absorption` array.
        This behaviour may change in the future.

        Can use `method` variations by providing the `var` parameter.


        Parameters
        ----------
        method : str
            Names or first letters of desired method.
        var : str, optional
            Name of the method variation, see `johnson_champoux`. The default is 'default'.

        Raises
        ------
        ValueError
            If some of the `method`'s required parameter is None
            or an unknown `method` is specified.

        Returns
        -------
        None.

        """
        if method.upper() in ['DB', 'DELANY-BAZLEY']:
            if not all([self.flowResistivity]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = delany_bazley(self.flowResistivity, self.air.density,
                                   self.air.soundSpeed, self.frequencies, var)

        elif method.upper() in ['R', 'RAY', 'RAYLEIGH']:
            if not all([self.flowResistivity, self.porosity]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = rayleigh(self.flowResistivity, self.air.density,
                              self.air.soundSpeed, self.porosity,
                              self.frequencies)

        elif method.upper() in ['BA', 'BIOT-ALLARD']:
            if not all([self.flowResistivity, self.porosity,
                        self.tortuosity, self.poreShape]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = biot_allard(self.flowResistivity, self.air.density, self.porosity,
                                 self.tortuosity, self.air.specHeatRatio, self.air.prandtl,
                                 self.air.atmPressure, self.poreShape, self.frequencies)

        elif method.upper() in ['JC', 'JOHNSON-CHAMPOUX']:
            if not all([self.flowResistivity, self.porosity, self.thermalLength,
                        self.tortuosity, self.viscousLength]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = johnson_champoux(self.flowResistivity, self.air.density,
                                      self.porosity, self.tortuosity,
                                      self.air.specHeatRatio, self.air.prandtl,
                                      self.air.atmPressure, self.viscousLength,
                                      self.thermalLength, self.air.viscosity,
                                      0 if not self.thermalPerm else self.thermalPerm,
                                      self.air.specHeatCP, self.frequencies, var)

        else:
            raise ValueError(f"Unknown method {method}.")
        self._zc = zc
        self._kc = kc
        self._absorp = absorption_coefficient(self.impedance, self.waveNum,
                                              self.thickness, self.air.impedance)
        return self.absorption
