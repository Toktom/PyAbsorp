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


def biot_allard(flow_resis, air_dens, poros, tortu, expans, prandtl,
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
            expans: int | float
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

    rhoEA = air_dens * tortu
    rhoEB = (flow_resis * poros) / (1j * omega * air_dens * tortu)
    rhoEC = (s * (-1j) ** 0.5) / 4
    rhoED = 2 / (s * (-1j) ** 0.5)
    rhoEE = ss.jv(1, s * (-1j) ** 0.5) / ss.jv(0, s * (-1j) ** 0.5)

    rhoE = rhoEA * (1 - rhoEB * ((rhoEC * rhoEE) / (1 - rhoED * rhoEE)))
    kEB = rhoEB / B
    kEC = rhoEC * B
    kED = rhoED / B
    kEE = ss.jv(1, s * B * (-1j) ** 0.5) / ss.jv(0, s * B * (-1j) ** 0.5)

    kE = (expans * atm) \
        / (expans - (expans - 1) / (1 - kEB * ((kEC * kEE) / (1 - kED * kEE))))

    # Charactheristic Impedance (zc) and the Wave Number (kc)
    zc = (kE * rhoE) ** 0.5
    kc = omega * (rhoE / kE) ** 0.5

    return zc, kc
