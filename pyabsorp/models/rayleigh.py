"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the rayleigh model.
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
    # Material Charactheristic Impedance (zc) and the Material Wave Number (kc)
    kc = (omega/sound_spd) * alpha
    zc = ((air_dens * sound_spd)/poros) * alpha
    return zc, kc
