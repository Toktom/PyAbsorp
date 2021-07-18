# -*- coding: utf-8 -*-
"""
Author: Michael Markus Ackermann
================================
Functions used to get the air properties.
"""
from typing import List


def air_properties(t: float, humidity: float, atmospheric_pressure: float,
                   kappa: float = 0.026, air_constant: float = 287.031,
                   water_constant: float = 461.521) -> List[float]:
    """Returns some of the air properties based on the ambient condition.

    Args:
        t (float): Air temperatue [°C].
        humidity (float): Air humidity [%].
        atmospheric_pressure (float): Atmospheric pressure [Pa].
        kappa (float, optional): W/(mK) air. Defaults to 0.026.
        air_constant (float, optional): Gas constant for air [J/K/kg].
                                        Defaults to 287.031.
        water_constant (float, optional): Gas constant for water vapor [J/K/kg].
                                Defaults to 461.521.

    Returns:
        c0 (float): The speed of the sound in the air [m/s].
        density (float): The air density [kg/(m^3)].
        impedance (float): Characteristic Impedance of the Air [Rayl].
        viscosity (float): The dynamic vicosity of the air (neta)[Ns/(m^2)].
        gama (float): Specific heat ratio [no units].
        prandtl_number (float):  Prandtl's number [no units].
        Cp (float): Spectfic Heat Constant Pressure [J/kg*K].
    """
    temp = ceilsius_to_kelvin(t)
    visc, p = viscosity(temp), pierce(temp)
    Cp = specific_heat_constant_pressure(temp)
    Cv = specific_heat_constant_volume(Cp)
    gama = specific_heat_ratio(Cp, Cv)
    prandtl_number = prandtl(visc, Cp, KAPPA=kappa)
    density = air_density(temp, humidity, atmospheric_pressure, p,
                          AIR_CONST=air_constant, WATER_CONST=water_constant)
    c0 = sound_speed(gama, atmospheric_pressure, density)
    impedance = air_impedance(c0, density)
    return c0, density, impedance, visc, gama, prandtl_number, Cp


def ceilsius_to_kelvin(t: float) -> float:
    """Converts the temperature from Ceilsius to Kelvin

    Args:
        t (float): Air temperatue [°C].

    Returns:
        float: Air temperature [K].
    """
    t += 273.16
    return t


def viscosity(t: float) -> float:
    """Calculates the dynamic vicosity of the air (neta (greek letter)).

    Args:
        t (float): Air temperature [K].

    Returns:
        float: Air dynamic vicosity [Ns/(m^2)].
    """
    return 7.72488e-8 * t - 5.95238e-11 * t**2 + 2.71368e-14 * t**3


def pierce(t: float) -> float:
    """Calculates Pierce.

    Args:
        t (float): Air temperature [K].

    Returns:
        float: Pierce.
    """
    return 0.0658 * t**3 - 53.7558 * t**2 + 14703.8127 * t - 1345485.0465


def specific_heat_constant_pressure(t: float) -> float:
    """Calculates the specific heat constant pressure.

    Args:
        t (float): Air temperature [K].

    Returns:
        float: Spectfic Heat Constant Pressure [J/kg*K].
    """
    return 4168.8 * (0.249679 - 7.55179e-5 * t + 1.69194e-7 * t**2 -
                     6.46128e-11 * t**3)


def specific_heat_constant_volume(Cp: float, AIR_CONST=287.031) -> float:
    """Calculates the specific heat constant volume for 260 K < T < 600 K

    Args:
        Cp (float):  Spectfic Heat Constant Pressure [J/kg*K].
        AIR_CONST (float, optional): Gas constant for air [J/K/kg].
            Defaults to 287.031.

    Returns:
        float: Spectfic Heat Constant Volume [J/kg/K].
    """
    return Cp - AIR_CONST


def prandtl(viscosity: float, Cp: float, KAPPA=0.026) -> float:
    """Calculatres the Prandtl number.

    Args:
        viscosity (float): Air dynamic vicosity [Ns/(m^2)].
        Cp (float): Spectfic Heat Constant Pressure [J/kg*K].
        KAPPA (float, optional): W/(mK) air. Defaults to 0.026.

    Returns:
        float: Prandtl number.
    """
    return viscosity * Cp / KAPPA


def specific_heat_ratio(Cp: float, Cv: float) -> float:
    """[summary]

    Args:
        Cp (float): Spectfic Heat Constant Pressure [J/kg*K].
        Cv (float): Spectfic Heat Constant Volume [J/kg/K].

    Returns:
        float: Specific heat ratio [no units].
    """
    return Cp / Cv


def air_density(t: float, humidity: float, atmospheric_pressure: float,
                pierce: float, AIR_CONST=287.031,
                WATER_CONST=461.521) -> float:
    """Calculates the air density.

    Args:
        t (float): Air temperature [K].
        humidity (float): Air humidity [%].
        atmospheric_pressure (float): Atmospheric pressure [Pa].
        pierce (float): Pierce.
        AIR_CONST (float, optional): Gas constant for air [J/K/kg].
                                        Defaults to 287.031.
        WATER_CONST (float, optional): Gas constant for water vapor [J/K/kg].
                                        Defaults to 461.521.

    Returns:
        float: Air density [kg/(m^3)]
    """
    return atmospheric_pressure / (AIR_CONST * t) \
        - (1 / AIR_CONST - 1 / WATER_CONST) * humidity / 100 * pierce / t


def sound_speed(gama: float, atmospheric_pressure: float,
                air_density: float) -> float:
    """Calculates the speed of the sound in air.

    Args:
        gama (float): Specific heat ratio [no units].
        atmospheric_pressure (float): Atmospheric pressure [Pa].
        air_density (float): Air density [kg/(m^3)]

    Returns:
        float: Sound speed in the air [m/s].
    """
    return (gama * atmospheric_pressure / air_density) ** 0.5


def air_impedance(sound_speed: float, air_density: float) -> float:
    """Calculates the Characteristic Impedance of the Air.

    Args:
        sound_speed (float): Sound speed in the air [m/s].
        air_density (float): Air density [kg/(m^3)]

    Returns:
        float: Characteristic Impedance of the Air [Rayl]
    """
    return sound_speed * air_density
