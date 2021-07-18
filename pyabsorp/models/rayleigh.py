# -*- coding: utf-8 -*-
"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the rayleigh model.
"""
import numpy as np


def rayleigh(flow_resis: float, air_dens: float, sound_spd: float,
             poros: float, freq=np.arange(100, 10001, 1)):
    """
    Returns through the Rayleigh Model the Material Charactheristic Impedance
    and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        sound_spd (float): The speed of sound in the air [m/s].
        poros (float): Material open porosity, between 0 and 1.
        freq (np.ndarray, optional): Array of frequencies.
            Defaults to np.arange(100, 10001, 1).

    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
    """
    omega = 2 * np.pi * freq
    alpha = (1 - (1j * poros * flow_resis) / (air_dens * omega)) ** 0.5
    # Material Charactheristic Impedance (zc) and the Material Wave Number (kc)
    kc = (omega/sound_spd) * alpha
    zc = ((air_dens * sound_spd)/poros) * alpha
    return zc, kc
