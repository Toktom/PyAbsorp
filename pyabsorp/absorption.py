"""
Author: Michael Markus Ackermann
================================
"""

from typing import List
import numpy as np


def absorption_coefficient(zc: List[complex], kc: List[complex], thickness: float,
                           z_air: float):
    """
    Returns the Sound Absorption Coefficient.
    NOTE 1: The values for 'zc' and 'kc' are already divided by the porosity.
    NOTE 2: This function only considers the normal incidence angle.

    Args:
        zc (List[complex]): Material Charactheristic Impedance.
        kc (List[complex]): Material Wave Number.
        thickness (float):  Material Thickness.
        z_air (float):  Air Characteristic Impedance.

    Returns:
        absorption (np. ndarray): Sound Absorption Coefficient [no units].
    """
    zs = -1j * (zc / np.tan(kc * thickness))  # Surface impedance (zs)
    vp = (zs - z_air) / (zs + z_air)  # Reflection coefficient (vp)
    return 1 - np.abs(vp) ** 2  # Sound Absorption Coefficient
