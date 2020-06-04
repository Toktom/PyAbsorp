"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the biot-allard model.
"""

import numpy as np
from scipy import special as ss

formats = {"circle": 1,
           "square": 1.07,
           "equi-tri": 1.11,
           "retang": 0.81}


def shear_wave(omega, flow_resis, poros, tortu, shape, air_dens):
    """
    Returns the Shear Wave number.

        Parameters:
        ----------
            omega: int | float | complex
                Angular frequency
            flow_ resis: int
                Resistivity of the material
            poros: float
                Porosity of the material
            tortu: float
                Tortuosity of the material
            shape: string
                Form factor for simple pores
            air_dens: int | float
                The air density

        Returns:
        --------
            s: int | float
                Shear wave number

    """
    c1 = formats[shape]
    num = 8 * omega * air_dens * tortu
    den = flow_resis * poros
    s = c1 * (num / den) ** 0.5

    return s


def biot_allard(flow_resis, air_dens, poros, tortu, gama, prandtl,
                atm, shape, freq=np.arange(100, 10001, 1)):
    """
    Returns through the Biot-Allard Model the Material Charactheristic
    Impedance and the Material Wave Number.

        Parameters:
        ----------
            flow_resis : int
                Resistivity of the material
            air_dens : int | float
                The air density
            poros : float
                Porosity of the material
            tortu: float
                Tortuosity of the material
            gama: int | float
                Ratio of specific heat
            prandtl: int | float
                Prandtl's number
            atm: int
                Atmospheric pressure
            shape: string
                Form factor for simple pores
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
    B = prandtl ** 0.5
    s = shear_wave(omega, flow_resis, poros, tortu, shape, air_dens)

    rho_part_a = air_dens * tortu
    rho_part_b = (flow_resis * poros) / (1j * omega * air_dens * tortu)
    rho_part_c = (s * (-1j) ** 0.5) / 4
    rho_part_d = 2 / (s * (-1j) ** 0.5)
    rho_part_e = ss.jv(1, s * (-1j) ** 0.5) / ss.jv(0, s * (-1j) ** 0.5)

    rho_ef = rho_part_a * (1 - rho_part_b * ((rho_part_c * rho_part_e) /
                           (1 - rho_part_d * rho_part_e)))

    k_part_a = rho_part_b / B
    k_part_b = rho_part_c * B
    k_part_c = rho_part_d / B
    k_part_d = ss.jv(1, s * B * (-1j) ** 0.5) / ss.jv(0, s * B * (-1j) ** 0.5)

    k_ef = (gama * atm) / (gama - (gama - 1) /
                           (1 - k_part_a * ((k_part_b * k_part_d) /
                            (1 - k_part_c * k_part_d))))

    # Changing from efficient to equivalent
    rho_eq = rho_ef / poros
    k_eq = k_ef / poros

    # Charactheristic Impedance (zc) and the Wave Number (kc)
    zc = (k_eq * rho_eq) ** 0.5
    kc = omega * (rho_eq / k_eq) ** 0.5

    return zc, kc
