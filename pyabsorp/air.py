"""
Author: Michael Markus Ackermann
================================
The function above describes the main purpouse of this file.
"""


def air_properties(temp, hum, atm):
    """
    Returns some of the air properties based on the ambient condition.

    Parameters:
    -----------
        temp: int | float
            Air temperatue (Celsius)
        hum : int | float
            Air humidity (%)
        atm : int | float
            Atmospheric pressure (Pa)

    Returns:
    --------
        sound_spd: int | float
            The speed of the sound in the air (m/s)
        air_dens: int | float
            The air density (kg/m³)
        z_air: int | float
            Characteristic Impedance of the Air
        viscos: int | float
            The dynamic vicosity of the air (a.k.a. neta (greek letter))(Ns/m²)
        gama: int | float
            Specific heat ratio (no units)
        prandtl: int | float
            Prandtl's number (fewly varies at typical air conditions)(no units)
        Cp: int | float
           Constant Pressure Spectfic Heat (J/kg·K)
    """
    KAPPA = 0.026  # W/(mK) air
    AIR_CONST = 287.031  # Gas constant for air (J/K/kg)
    WATER_CONST = 461.521  # Gas constant for water vapor (J/K/kg)
    temp += 273.16  # Converts the temperature from Ceilsius to Kelvin
    pierce = 0.0658 * temp**3 - 53.7558 * temp**2 \
        + 14703.8127 * temp - 1345485.0465  # Pierce
    viscos = 7.72488e-8 * temp - 5.95238e-11 * temp**2 \
        + 2.71368e-14 * temp**3  # Absolute (or dynamic) viscosity
    Cp = 4168.8 * (0.249679 - 7.55179e-5 * temp
                   + 1.69194e-7 * temp**2 - 6.46128e-11 * temp**3)
    # Constant Pressure Spectfic Heat
    Cv = Cp - AIR_CONST  # Constant Volume Specific Heat (J/kg/K) for 260 K < T < 600 K
    prandtl = viscos * Cp / KAPPA  # Prandtl number
    gama = Cp / Cv  # Specific heat ratio
    air_dens = atm / (AIR_CONST * temp) - (1 / AIR_CONST - 1 / WATER_CONST) * \
        hum / 100 * pierce / temp  # Air density
    sound_spd = (gama * atm / air_dens) ** 0.5  # Speed of the sound
    z_air = sound_spd * air_dens  # Characteristic Impedance of the Air
    return sound_spd, air_dens, z_air, viscos, gama, prandtl, Cp


class AirProperties(object):
    """Air acoustical properties object interface."""

    def __init__(self, temp: float or int = 20, hum: float or int = 50,
                 atm: float or int = 101325):
        """
        Air properties for acoustical parameters.

        Parameters
        ----------
        temp : float or int, optional
            Environment temperature, in Celsius degrees [°C]. The default is 20.
        hum : float or int, optional
            Relative humidity, as percentage [%]. The default is 50.
        atm : float or int, optional
            Atmospheric pressure, in Pascals [Pa]. The default is 101325.

        Returns
        -------
        None.

        """
        self.temperature = temp
        self.humidity = hum
        self.atmPressure = atm
        self.calculate_properties()
        return

    def calculate_properties(self):
        """
        Calculate the air acoustically relevant properties.

        The `AirProperties` basic values, `temperature`, `humidity` and
        `atmPressure` can be set at will, but this method must be explicitly
        called to update the properties values.

        Returns
        -------
        None.

        """
        props = air_properties(self.temperature, self.humidity, self.atmPressure)
        self._soundSpeed = props[0]
        self._density = props[1]
        self._impedance = props[2]
        self._viscosity = props[3]
        self._specHeatRatio = props[4]
        self._prandtl = props[5]
        self._specHeatCP = props[6]
        self._specHeatCV = self.specHeatCP / self.specHeatRatio
        return

    @property
    def temperature(self):
        """Air temperature in degree Celsius [°C]."""
        return self._temp

    @temperature.setter
    def temperature(self, temp: float or int):
        if type(temp) not in [float, int]:
            raise TypeError("Temperature must be a number.")
        elif temp < 0 or temp > 50:
            raise ValueError("Temperature must be between 0 and 50 °C.")
        self._temp = temp
        return

    @property
    def humidity(self):
        """Air relative humidity in percentage [%]."""
        return self._hum

    @humidity.setter
    def humidity(self, hum):
        if type(hum) not in [float, int]:
            raise TypeError("Humidity must be a number.")
        elif hum < 0 or hum > 50:
            raise ValueError("Humidity must be between 0 and 100 %.")
        self._hum = hum
        return

    @property
    def atmPressure(self):
        """Atmospherical pressure in Pascals [Pa]."""
        return self._atm

    @atmPressure.setter
    def atmPressure(self, atm):
        if type(atm) not in [float, int]:
            raise TypeError("Atmospheric pressure must be a number.")
        elif atm < 90e3 or atm > 115e3:
            raise ValueError("Atmospheric pressure must be between 90 and 115 kPa.")
        self._atm = atm
        return

    @property
    def soundSpeed(self):
        """Sound speed in meters per second [m/s]."""
        return self._soundSpeed

    @property
    def density(self):
        """Specific, or volumetric, density in kilograms per cubic meter [kg/m³]."""
        return self._density

    @property
    def impedance(self):
        """Characteristic impedance in rayls per square meter [rayl/m²]."""
        return self._impedance

    @property
    def specHeatCP(self):
        """Specific heat at Constant Pressure [J/kgK]."""
        return self._specHeatCP

    @property
    def specHeatCV(self):
        """Specific heat at Constant Volume [J/kgK]."""
        return self._specHeatCV

    @property
    def specHeatRatio(self):
        """
        Ratio between specific heat at constant pressure
        and specific heat at constant volume [-].
        """
        return self._specHeatRatio

    @property
    def viscosity(self):
        """Dynamic viscosity [Ns/m²]."""
        return self._viscosity

    @property
    def prandtl(self):
        """Prandtl number [-]."""
        return self._prandtl
