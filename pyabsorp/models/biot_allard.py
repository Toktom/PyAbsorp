"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the biot-allard model.
"""
import numpy as np
from typing import List
from scipy import special as ss


formats = {"circle": 1,
           "square": 1.07,
           "equilateral triangular": 1.11,
           "retangular": 0.81}


def shear_wave(omega: List[float], flow_resis: float, poros: float, tortu: float,
               shape: str, air_dens: float) -> float:
    """
    Returns the Shear Wave number.

    Args:
        omega (List[float]): Array of angular frequencies.
        flow_resis (float): Static flow resistivity of the material  [Ns/(m^4)].
        poros (float): Material open porosity, between 0 and 1.
        tortu (float): Material tortuosity.
        shape (str): Form factor for simple pores, must be a 'circle',
            'square', 'equilateral triangular' or 'retangular'.
        air_dens (float): Air density [kg/(m^3)].

    Returns:
        float: Shear Wave number.
    """
    c1 = formats[shape]
    num = 8 * omega * air_dens * tortu
    den = flow_resis * poros
    return c1 * (num / den) ** 0.5


def biot_allard(flow_resis: float, air_dens: float, poros: float, tortu: float,
                gama: float, prandtl: float, atm: float, shape: str,
                freq=np.arange(100, 10001, 1)):
    """
    Returns through the Biot-Allard Model the Material Charactheristic
    Impedance and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material  [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        poros (float): Material open porosity, between 0 and 1.
        tortu (float): Material tortuosity.
        gama (float): Specific heat ratio [no units].
        prandtl (float): Prandtl's number.
        atm (float): Atmospheric pressure [Pa].
        shape (str) : Form factor for simple pores, must be a 'circle',
            'square', 'equilateral triangular' or 'retangular'.
        freq (np.ndarray, optional): Array of frequencies.
            Defaults to np.arange(100, 10001, 1).


    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
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
