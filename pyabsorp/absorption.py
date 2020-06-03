"""
Author: Michael Markus Ackermann
================================
"""

import numpy as np


def absorption_coefficient(zc, kc, d, z_air):
    """
    Returns the Sound Absorption Coefficient.
    NOTE: This function only considers the normal incidence angle.

        Parameters:
        -----------
            zc : int | float | complex
                Material Charactheristic Impedance
            kc : int | float | complex
                Material Wave Number
            d : float
                Material Thickness
            z_air : int | float
                Air Characteristic Impedance
            poros: float
                Porosity of the Material

        Returns:
        --------
            absorption : int | ndarray
                Sound Absorption Coefficient of the Material
    """
    zs = -1j * (zc / np.tan(kc * d))  # Surface impedance (zs)
    vp = (zs - z_air) / (zs + z_air)  # Reflection coefficient (vp)
    absorption = 1 - np.abs(vp) ** 2  # Sound Absorption Coefficient
    return absorption
