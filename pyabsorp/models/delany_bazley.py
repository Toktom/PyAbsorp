"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the delany-bazley model.
"""

import numpy as np


def delany_bazley(flow_resis, air_dens, sound_spd,
                  freq=np.arange(100, 10001, 1), var='default'):
    """
    Returns through the Delany-Bazley Model the Material Charactheristic
    Impedance and the Material Wave Number.

        Parameters:
        ----------
            flow_resis : int
                Resistivity of the material
            air_dens : int | float
                The air density
            sound_spd : int | float
                The speed of the sound
            freq : ndarray
                A range of frequencies
                NOTE: default range goes from 100 [Hz] to 10 [kHz].
            var : string
                The variation of the Delany-Bazley Model
                Variations availabe:
                -'default'           -> Delany-Bazley
                -'miki'              -> Delany-Bazley-Miki
                -'allard-champoux'   -> Delany-Bazley-Allard-Champoux

        Returns:
        -------
            zc : int | float | complex
                Material Charactheristic Impedance
            kc : int | float | complex
                Material Wave Number
    """
    if var == 'default':  # Original Delany Bazley
        R = 1 + 9.08 * ((1e3 * freq / flow_resis) ** -0.75)
        X = -11.9 * ((1e3 * freq / flow_resis) ** -0.73)
        # Charactheristic Impedance (zc)
        zc = sound_spd * air_dens * (R + 1j * X)

        alpha = 1 + 10.8 * (1e3 * freq / flow_resis) ** -0.7
        beta = -10.3 * (1e3 * freq / flow_resis) ** -0.59
        # Wave Number (kc)
        kc = (2 * np.pi * freq / sound_spd) * (alpha + 1j * beta)

    elif var == 'miki':  # Miki variation
        R = 1 + 5.50 * ((1e3 * freq / flow_resis) ** -0.632)
        X = -8.43 * ((1e3 * freq / flow_resis) ** -0.632)
        # Charactheristic Impedance (zc)
        zc = sound_spd * air_dens * (R + 1j * X)

        alpha = 1 + 7.81 * (1e3 * freq / flow_resis) ** -0.618
        beta = -11.41 * (1e3 * freq / flow_resis) ** -0.618
        # Wave Number (kc)
        kc = (2 * np.pi * freq / sound_spd) * (alpha + 1j * beta)

    elif var == 'allard-champoux':  # Allard and Champoux variation
        R = 1 + 0.0571 * (((air_dens*freq) / flow_resis) ** -0.754)
        X = -0.0870 * (((air_dens*freq) / flow_resis) ** -0.732)
        # Charactheristic Impedance (zc)
        zc = sound_spd * air_dens * (R + 1j * X)

        alpha = 1 + 0.0978 * ((air_dens*freq) / flow_resis) ** -0.7
        beta = -0.1890 * ((air_dens*freq) / flow_resis) ** -0.595
        # Wave Number (kc)
        kc = (2 * np.pi * freq / sound_spd) * (alpha + 1j * beta)

    return zc, kc
