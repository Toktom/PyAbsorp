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
