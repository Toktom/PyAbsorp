# -*- coding: utf-8 -*-
"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the delany-bazley model.
"""
import numpy as np
from typing import Tuple


def delany_bazley(flow_resis: float, air_dens: float, sound_spd: float,
                  freq: np.ndarray = np.arange(100, 10001, 1),
                  var: str = 'default') -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns through the Delany-Bazley Model the Material Charactheristic
    Impedance and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        sound_spd (float): The speed of sound in the air [m/s].
        freq (np.ndarray, optional): Array of frequencies.
            Defaults to np.arange(100, 10001, 1).
        var (str, optional): Model variation. Defaults to 'default'.
        Model variations availabe:
            -'default'          -> Original Delany Bazley
            -'miki'             -> Miki variation
            -'allard-champoux'  -> Allard and Champoux variation


    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
    """
    if var == 'allard-champoux':
        (zc, kc) = __allard_champoux(flow_resis, air_dens, sound_spd, freq)

    elif var == 'default':
        (zc, kc) = __default(flow_resis, air_dens, sound_spd, freq)

    elif var == 'miki':
        (zc, kc) = __miki(flow_resis, air_dens, sound_spd, freq)
    else:
        raise ValueError(f"'var = {var}' isn't a valid value.")
    return (zc, kc)


def __allard_champoux(flow_resis: float, air_dens: float, sound_spd: float,
                      freq: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns through the Delany-Bazley-Allard-Champoux Model the Material
    Charactheristic Impedance and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        sound_spd (float): The speed of sound in the air [m/s].
        freq (np.ndarray, optional): Array of frequencies.
            Defaults to np.arange(100, 10001, 1).

    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
    """
    # BIG NOTE: R stands for REAL part and X for IMAGINARY part
    R = 1 + 0.0571 * (((air_dens*freq) / flow_resis) ** -0.754)
    X = -0.0870 * (((air_dens*freq) / flow_resis) ** -0.732)
    # Charactheristic Impedance (zc)
    zc = sound_spd * air_dens * (R + 1j * X)

    alpha = 1 + 0.0978 * ((air_dens*freq) / flow_resis) ** -0.7
    beta = -0.1890 * ((air_dens*freq) / flow_resis) ** -0.595
    # Wave Number (kc)
    kc = (2 * np.pi * freq / sound_spd) * (alpha + 1j * beta)
    return (zc, kc)


def __default(flow_resis: float, air_dens: float, sound_spd: float,
              freq: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns through the Delany-Bazley Model the Material
    Charactheristic Impedance and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        sound_spd (float): The speed of sound in the air [m/s].
        freq (np.ndarray, optional): Array of frequencies.
            Defaults to np.arange(100, 10001, 1).

    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
    """
    # BIG NOTE: R stands for REAL part and X for IMAGINARY part
    R = 1 + 9.08 * ((1e3 * freq / flow_resis) ** -0.75)
    X = -11.9 * ((1e3 * freq / flow_resis) ** -0.73)
    # Charactheristic Impedance (zc)
    zc = sound_spd * air_dens * (R + 1j * X)

    alpha = 1 + 10.8 * (1e3 * freq / flow_resis) ** -0.7
    beta = -10.3 * (1e3 * freq / flow_resis) ** -0.59
    # Wave Number (kc)
    kc = (2 * np.pi * freq / sound_spd) * (alpha + 1j * beta)
    return (zc, kc)


def __miki(flow_resis: float, air_dens: float, sound_spd: float,
           freq: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns through the Delany-Bazley-Miki Model the Material
    Charactheristic Impedance and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        sound_spd (float): The speed of sound in the air [m/s].
        freq (np.ndarray, optional): Array of frequencies.
            Defaults to np.arange(100, 10001, 1).

    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
    """
    # BIG NOTE: R stands for REAL part and X for IMAGINARY part
    R = 1 + 5.50 * ((1e3 * freq / flow_resis) ** -0.632)
    X = -8.43 * ((1e3 * freq / flow_resis) ** -0.632)
    # Charactheristic Impedance (zc)
    zc = sound_spd * air_dens * (R + 1j * X)

    alpha = 1 + 7.81 * (1e3 * freq / flow_resis) ** -0.618
    beta = -11.41 * (1e3 * freq / flow_resis) ** -0.618
    # Wave Number (kc)
    kc = (2 * np.pi * freq / sound_spd) * (alpha + 1j * beta)
    return (zc, kc)
