"""
Author: Michael Markus Ackermann
================================
Air class.
"""
from dataclasses import field, dataclass
from pyabsorp.air import air_properties


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
        water_constant (float, optional): Gas constant for water vapor (J/K/kg).
                                Defaults to 461.521.

    Returns:
        None.

    """
    temperature: float = field(default=20, metadata=dict(
        title="Enverioment temperature", description="Should be in Celsius degrees [°C]"))

    humidity: float = field(default=50, metadata=dict(
        title="Relative humidity", description="Should be in as percentage [%]"))

    atmospheric_pressure: float = field(default=101325, metadata=dict(
        title="Atmospheric pressure", description="Should be in Pascals [Pa]"))

    kappa: float = field(default=0.026, metadata=dict(
        title="Kappa", description="Constant"))

    air_constant: float = field(default=287.031, metadata=dict(
        title="Gas constant for air", description="Should be in [J/K/kg]"))

    water_constant: float = field(default=461.521, metadata=dict(
        title="Gas constant for water steam", description="Should be in [J/K/kg]"))

    def __post_init__(self):
        properties = air_properties(self.temperature, self.humidity, self.atmospheric_pressure,
                                    self.kappa, self.air_constant, self.water_constant)
        self.speed = properties[0]
        self.density = properties[1]
        self.impedance = properties[2]
        self.viscosity = properties[3]
        self.specific_heat_ratio = properties[4]
        self.prandtl = properties[5]
        self.specific_heat_cp = properties[6]


    def __repr__(self) -> str:
        r = f'Air(temperature={self.temperature}, humidity={self.humidity}, atmospheric_pressure={self.atmospheric_pressure})'
        return r
