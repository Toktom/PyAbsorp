"""

=================
@Author: Michael Markus Ackermann
"""


def air_properties(temp, hum, atm):
    """
    Returns some of the air properties based on the ambient condition.

    Parameters:
    -----------
        temp: int | float
            Air temperatue
        hum : int | float
            Air humidity
        atm : int | float
            Atmospheric pressure

    Returns:
    --------
        sound_spd: int | float
            The speed of the sound in the air
        air_dens: int | float
            The air density
        c_imp_air: int | float
            Characteristic Impedance of the Air
        viscos: int | float
            The dynamic vicosity of the air (a.k.a. neta (greek letter))
        expans: int | float
            Ratio of specific heat
        prandtl: int | float
            Prandtl's number
        Cp: int | float
            Specific heat capacity
    """
    KAPPA = 0.026  # W/(mK) air
    AIR_CONST = 287.031  # J/K/kg
    WATER_CONST = 461.521  # J/K/kg
    temp += 273.16  # K#
    pierce = 0.0658 * temp**3 - 53.7558 * temp**2 \
        + 14703.8127 * temp - 1345485.0465
    viscos = 7.72488e-8 * temp - 5.95238e-11 * temp**2 \
        + 2.71368e-14 * temp**3
    Cp = 4168.8 * (0.249679 - 7.55179e-5 * temp
                   + 1.69194e-7 * temp**2 - 6.46128e-11 * temp**3)
    # Volume Specific Heat Constant (J/kg/K) for 260 K < T < 600 K
    vol_Const = Cp - AIR_CONST
    prandtl = viscos * Cp / KAPPA  # prandtl number
    expans = Cp / vol_Const  # k
    air_dens = atm / (AIR_CONST * temp) - (1 / AIR_CONST - 1 / WATER_CONST) * \
        hum / 100 * pierce / temp  # Air density
    sound_spd = (expans * atm / air_dens) ** 0.5  # Speed of the sound
    c_imp_air = sound_spd * air_dens  # Characteristic Impedance of the Air
    return sound_spd, air_dens, c_imp_air, viscos, expans, prandtl, Cp
