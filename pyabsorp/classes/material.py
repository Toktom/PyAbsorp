# -*- coding: utf-8 -*-
"""
Author: João Vitor Gutkoski Paes
Coauthor: Michael Markus Ackermann
================================
Provide basic interface to handle a single material being studied.
"""
import numpy as np
from dataclasses import field, dataclass
from pyabsorp.classes import Air
from pyabsorp.functions import absorption_coefficient
from pyabsorp.models import delany_bazley, rayleigh, biot_allard, \
    johnson_champoux
import pyabsorp.classes.__classes_metadata as cm


@dataclass
class Material:
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

    Args:

    air (Air): Air object from pyabsorp.Air.
    thickness (float): Material thickness [m].
    porosity (float) [LP]: Material open porosity, between 0 and 1.
    tortuosity (float) [LP]: Material tortuosity.
    flow_resistivity (float) [LP]: Satic flow resistivity [Ns/(m^4)].
    thermal_length (float) [LP]: Thermal characteristic length [m].
    viscosity_length (float) [LP]: Viscous characteristic length [m].
    thermal_permeability (float) [LP]: Static thermal permeability [m^2].
    shape (str) [LP]: The shape of the pore, must be a 'circle', 'square',
        'equilateral triangular' or 'retangular'.
    frequencies (np.ndarray): Array of frequencies.

    Returns:
        None.
    """

    air: Air = field(default=None, metadata=dict(
        title=cm.air_md[0], description=cm.air_md[1]))
    thickness: float = field(default=None, metadata=dict(
        title=cm.thick_md[0], description=cm.thick_md[1]))
    porosity: float = field(default=None, metadata=dict(
        title=cm.poro_md[0], description=cm.poro_md[1]))
    tortuosity: float = field(default=None, metadata=dict(
        title=cm.tortu_md[0], description=cm.tortu_md[1]))
    flow_resistivity: float = field(default=None, metadata=dict(
        title=cm.flow_md[0], description=cm.flow_md[1]))
    thermal_length: float = field(default=None, metadata=dict(
        title=cm.th_l_md[0], description=cm.th_l_md[1]))
    viscosity_length: float = field(default=None, metadata=dict(
        title=cm.vs_l_md[0], description=cm.vs_l_md[1]))
    thermal_permeability: float = field(default=None, metadata=dict(
        title=cm.th_p_md[0], description=cm.th_p_md[1]))
    shape: str = field(default=None, metadata=dict(
        title=cm.shape_md[0], description=cm.shape_md[1]))
    frequencies: np.ndarray = field(default=np.arange(100, 10001, 1),
                                    metadata=dict(title=cm.freq_md[0],
                                    description=cm.freq_md[1]))

    __annotations__ = {
        'air': Air,
        'thickness': float,
        'porosity': float,
        'tortuosity': float,
        'flow_resistivity': float,
        'thermal_length': float,
        'viscosity_length': float,
        'thermal_permeability': float,
        'shape': str,
        'frequencies': np.ndarray
    }

    def estimate_absorption(self, frequencies: np.ndarray, method: str,
                            var: str = 'default'):
        """
        Estimate material absorption based on `method`.

        The material will hold the resulting `absorption` coefficients and the
        respective characteristic `impedance` and wave number (`waveNum`).

        Only the result of one call can be held. This means that a comparison
        between methods, or method variations must save separatedly each
        `absorption` array. This behaviour may change in the future.

        Can use `method` variations by providing the `var` parameter.


        Args:
            frequencies (np.ndarray): Array of frequencies used to estimate
                `impedance`, `waveNum` and `absorption`.
            method (str): Names or first letters of desired method, e.g.
                'rayleigh' for Rayleigh, or 'jc' for Johnson-Champoux.
            var (str, optional): Name of the method variation, see
                `johnson_champoux`. The default is 'default'.

        Raises:
            ValueError: If some of the `method`'s required parameter is None
                or an unknown `method` is specified.

        Returns:
            None.

        """
        if method.upper() in ['DB', 'DELANY-BAZLEY']:
            if not all([self.flow_resistivity]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = delany_bazley(self.flow_resistivity, self.air.density,
                                   self.air.speed, frequencies, var)

        elif method.upper() in ['R', 'RAY', 'RAYLEIGH']:
            if not all([self.flow_resistivity, self.porosity]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = rayleigh(self.flow_resistivity, self.air.density,
                              self.air.speed, self.porosity,
                              frequencies)

        elif method.upper() in ['BA', 'BIOT-ALLARD']:
            if not all([self.flow_resistivity, self.porosity,
                        self.tortuosity, self.shape]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = biot_allard(self.flow_resistivity, self.air.density,
                                 self.porosity, self.tortuosity,
                                 self.air.specific_heat_ratio,
                                 self.air.prandtl,
                                 self.air.atmospheric_pressure,
                                 self.shape, frequencies)

        elif method.upper() in ['JC', 'JOHNSON-CHAMPOUX']:
            if not all([self.flow_resistivity, self.porosity,
                        self.thermal_length, self.tortuosity,
                        self.viscosity_length]):
                raise ValueError("Some material parameters are not defined.")

            zc, kc = johnson_champoux(self.flow_resistivity, self.air.density,
                                      self.porosity, self.tortuosity,
                                      self.air.specific_heat_ratio,
                                      self.air.prandtl,
                                      self.air.atmospheric_pressure,
                                      self.viscosity_length,
                                      self.thermal_length, self.air.viscosity,
                                      0 if not self.thermal_permeability else self.thermal_permeability,
                                      self.air.specific_heat_cp, frequencies,
                                      var)

        else:
            raise ValueError(f"Unknown method {method}.")
        self.frequencies = frequencies
        self.impedance = zc
        self.wave_number = kc
        self.absorption = absorption_coefficient(self.impedance,
                                                 self.wave_number,
                                                 self.thickness,
                                                 self.air.impedance)
        return self.absorption

    def __repr__(self) -> str:
        return f"""Material(\nthickness={self.thickness},
            \nporosity={self.porosity},\ntortuosity={self.tortuosity}
            \nflow_resistivity={self.flow_resistivity},\nthermal_length={self.thermal_length},
            \nviscosity_length={self.viscosity_length}
            \nthermal_permeability={self.thermal_permeability},\nshape={self.shape},
            \nfrequencies={self.frequencies})"""
