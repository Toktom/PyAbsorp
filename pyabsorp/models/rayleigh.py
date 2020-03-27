"""

=================
@Author: Michael Markus Ackermann
"""

import numpy as np


def rayleigh(flow_resis, air_dens, sound_spd,
             poros, freq=np.arange(100, 10001, 1)):
    """
    Returns through the Rayleigh Model the Material Charactheristic Impedance
    and the Material Wave Number.

        Parameters:
        ----------
            flow_resis : int
                Resistivity of the material
            air_dens : int | float
                The air density
            sound_spd : int | float
                The speed of the sound
            poros : float
                Porosity of the material
            freq : ndarray
                A range of frequencies
                NOTE: default range goes from 100 [Hz] to 10 [kHz].

        Returns:
        -------
            zc : int | float | complex
                Material Charactheristic Impedance
            kc : int | float | complex
                Material Wave Number
    """
    omega = 2 * np.pi * freq
    alpha = (1 - (1j * poros * flow_resis) / (air_dens * omega)) ** 0.5
    kc = (omega/sound_spd) * alpha
    zc = ((air_dens * sound_spd)/poros) * alpha
    return zc, kc


def rayleigh_absorption(zc, kc, d, c_imp_air):
    """
    Returns the Sound Absorption Coefficient for the Rayleigh Model.
    NOTE: Only use it with the rayleigh function.

        Parameters:
        -----------
            zc : int | float | complex
                Material Charactheristic Impedance
            kc : int | float | complex
                Material Wave Number
            d : float
                Material Thickness
            c_imp_air : int | float
                Air Characteristic Impedance

        Returns:
        --------
            absorption : int | ndarray
                Sound Absorption Coefficient of the Material
    """
    zs = -1j * (zc / np.tan(kc * d))
    reflex = (zs - c_imp_air) / (zs + c_imp_air)
    absorption = 1 - np.abs(reflex) ** 2
    return absorption
