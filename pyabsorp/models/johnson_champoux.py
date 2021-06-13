"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the johnson-champoux model.
"""
import numpy as np


def johnson_champoux(flow_resis: float, air_dens: float, poros: float, tortu: float,
                     gama: float, prandtl: float, atm: float, visc: float, therm: float,
                     neta: float, therm_perm: float = 0, Cp: float = 0,
                     freq: np.ndarray = np.arange(100, 10001, 1), var: str = 'default'):
    """Returns through the Johnson-Champoux Model the Material Charactheristic
    Impedance and the Material Wave Number.

    Args:
        flow_resis (float): Static flow resistivity of the material  [Ns/(m^4)].
        air_dens (float): Air density [kg/(m^3)].
        poros (float): Material open porosity, between 0 and 1.
        tortu (float): Material tortuosity.
        gama (float): Specific heat ratio [no units].
        prandtl (float): Prandtl's number [no units].
        atm (float): Atmospheric pressure [Pa].
        visc (float): Viscous characteristic length [m].
        therm (float): Thermal characteristic length [m].
        neta (float):  The dynamic vicosity of the air (neta (greek letter))[Ns/(m^2)].
        therm_perm (float, optional): [description]. Defaults to 0.
        Cp (float, optional): Spectfic Heat Constant Pressure [J/kg*K]. Defaults to 0.
        freq (np.ndarray, optional): Array of frequencies. 
            Defaults to np.arange(100, 10001, 1).
        var (str, optional): Model variation. Defaults to 'default'.
        Model variations availabe:
            -'default'           -> Johnson-Champoux
            -'allard'            -> Johnson-Champoux-Allard
            -'lafarge'           -> Johnson-Champoux-Allard-Lafarge


    Returns:
        zc (np.ndarray): Material Charactheristic Impedance.
        kc (np.ndarray): Material Wave Number.
    """
    omega = 2 * np.pi * freq

    if var == 'default':
        # rho_ef
        rho_ef_part_a = 4 * omega * air_dens * neta * (tortu**2)
        rho_ef_part_b = ((flow_resis**2) * (poros**2) * (visc**2))
        rho_ef_part_c = (1 + 1j * (rho_ef_part_a / rho_ef_part_b)) ** 0.5
        rho_ef_part_d = 1j * omega * air_dens * tortu
        rho_ef_part_e = 1 + ((poros*flow_resis)/rho_ef_part_d) * rho_ef_part_c

        rho_ef = air_dens * tortu * rho_ef_part_e

        # k_ef
        k_ef_part_a = omega * prandtl * air_dens * (therm ** 2)
        k_ef_part_b = (1 + 1j * (k_ef_part_a / (16 * neta))) ** 0.5
        k_ef_part_c = (1 - ((1j * 8 * neta)/k_ef_part_a) * k_ef_part_b) ** -1
        k_ef_part_d = gama - (gama - 1) * k_ef_part_c
        k_ef = (gama * atm) / k_ef_part_d

        # Changing from efficient to equivalent
        rho_eq = rho_ef / poros
        k_eq = k_ef / poros

        # Charactheristic Impedance (zc) and the Wave Number (kc)
        kc = omega * ((rho_eq / k_eq) ** 0.5)
        zc = (k_eq * rho_eq) ** 0.5

        return zc, kc

    if var == 'allard':
        # rho_ef
        rho_ef_part_a = 4 * air_dens * (tortu**2) * neta * omega
        rho_ef_part_b = (flow_resis**2) * (poros**2) * (visc**2)
        rho_ef_part_c = (1 + 1j * (rho_ef_part_a/rho_ef_part_b)) ** 0.5
        rho_ef_part_d = (flow_resis * poros) / (1j * air_dens * tortu * omega)
        rho_ef_part_e = 1 + rho_ef_part_d * rho_ef_part_c
        rho_ef = (air_dens * tortu) * rho_ef_part_e

        # k_ef
        k_ef_part_a = 4 * air_dens * (tortu**2) * neta * omega
        k_ef_part_b = (flow_resis**2) * (poros**2) * (therm**2)
        k_ef_part_c = (1 + 1j * (k_ef_part_a/k_ef_part_b)) ** 0.5
        k_ef_part_d = (flow_resis * poros)/(air_dens * prandtl * tortu * omega)
        k_ef_part_e = (1 - 1j * k_ef_part_d * k_ef_part_c) ** -1
        k_ef_part_f = gama - (gama - 1) * k_ef_part_e
        k_ef = (gama * atm) / k_ef_part_f

        # Changing from efficient to equivalent
        rho_eq = rho_ef / poros
        k_eq = k_ef / poros

        # Charactheristic Impedance (zc) and the Wave Number (kc)
        kc = omega * ((rho_eq / k_eq) ** 0.5)
        zc = (k_eq * rho_eq) ** 0.5

        return zc, kc

    if var == "lafarge":
        KAPPA = 0.026  # W/(mK) air
        # Static thermal permeability

        # Common part used to obtain rho_ef and k_ef
        common_part = 4 * omega * air_dens

        # rho_ef
        rho_ef_part_a = common_part * neta * (tortu ** 2)
        rho_ef_part_b = (flow_resis ** 2) * (visc ** 2) * (poros ** 2)
        rho_ef_part_c = (1 + 1j * (rho_ef_part_a/rho_ef_part_b)) ** 0.5
        rho_ef_part_d = (flow_resis * poros) / (1j*omega*air_dens*tortu)
        rho_ef_part_e = 1 + (rho_ef_part_d * rho_ef_part_c)
        rho_ef = (air_dens * tortu) * rho_ef_part_e

        # k_ef
        k_ef_part_a = common_part * Cp * (therm_perm ** 2)
        k_ef_part_b = KAPPA * (therm ** 2) * (poros ** 2)
        k_ef_part_c = (1 + 1j * (k_ef_part_a / k_ef_part_b)) ** 0.5
        k_ef_part_d = (poros * KAPPA) / (therm_perm * Cp * air_dens * omega)
        k_ef_part_e = (1 - 1j * (k_ef_part_d * k_ef_part_c)) ** -1
        k_ef_part_f = gama - (gama - 1) * k_ef_part_e
        k_ef = (gama * atm) / k_ef_part_f

        # Changing from efficient to equivalent
        rho_eq = rho_ef / poros
        k_eq = k_ef / poros

        # Charactheristic Impedance (zc) and the Wave Number (kc)
        kc = omega * ((rho_eq / k_eq) ** 0.5)
        zc = (k_eq * rho_eq) ** 0.5

        return zc, kc
