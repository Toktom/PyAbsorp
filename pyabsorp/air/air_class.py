# -*- coding: utf-8 -*-
"""
Author: Michael Markus Ackermann
================================
Air class.
"""
from dataclasses import field, dataclass
from pyabsorp.air import air_properties
import pyabsorp.utils.__classes_metadata as cm


@dataclass
class Air:
    """
        Air properties class used in acoustics.

    Args:
        temperature (float): Air temperatue (Celsius).
        humidity (float): Air humidity (%).
        atmospheric_pressure (float): Atmospheric pressure (Pa).
        kappa (float, optional): W/(mK) air. Defaults to 0.026.
        air_constant (float, optional): Gas constant for air (J/K/kg).
                                        Defaults to 287.031.
        water_constant (float, optional): Gas constant for water steam (J/K/kg).
                                Defaults to 461.521.

    Returns:
        None.

    """
    temperature: float = field(default=20, metadata=dict(
        title=cm.temp_md[0], description=cm.temp_md[1]))

    humidity: float = field(default=50, metadata=dict(
        title=cm.hum_md[0], description=cm.hum_md[1]))

    atmospheric_pressure: float = field(default=101325, metadata=dict(
        title=cm.atm_md[0], description=cm.atm_md[1]))

    kappa: float = field(default=0.026, metadata=dict(
        title=cm.kappa_md[0], description=cm.kappa_md[1]))

    air_constant: float = field(default=287.031, metadata=dict(
        title=cm.air_c_md[0], description=cm.air_c_md[1]))

    water_constant: float = field(default=461.521, metadata=dict(
        title=cm.water_c_md[0], description=cm.water_c_md[1]))

    def __post_init__(self):
        properties = air_properties(self.temperature, self.humidity,
                                    self.atmospheric_pressure, self.kappa,
                                    self.air_constant, self.water_constant)
        self.speed = properties[0]
        self.density = properties[1]
        self.impedance = properties[2]
        self.viscosity = properties[3]
        self.specific_heat_ratio = properties[4]
        self.prandtl = properties[5]
        self.specific_heat_cp = properties[6]

    def __repr__(self) -> str:
        return f'Air(temperature={self.temperature}, humidity={self.humidity},\
            atmospheric_pressure={self.atmospheric_pressure})'
