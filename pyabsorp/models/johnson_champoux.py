"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the johnson-champoux model.
"""


import numpy as np


def johnson_champoux(flow_resis, air_dens, poros, tortu, gama,
                     prandtl, atm, visc, therm, neta, therm_perm=0, Cp=0,
                     freq=np.arange(100, 10001, 1), var='default'):
    """
    Returns through the Johnson-Champoux Model the Material Charactheristic
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
            visc: int | float
                Viscous characteristic length
            therm: int | float
                Thermal characteristic length
            neta: int | float
                Dynamic vicosity of air
            therm_perm: int | float
                Static thermal permeability
            Cp: int | float
                Specific heat capacity at constant pressure
                NOTE: Only used int Johnson-Champoux-Allard model.
            freq : ndarray
                A range of frequencies
                NOTE: default range goes from 100 [Hz] to 10 [kHz].
            var : string
                The variation of the Johnson-Champoux Model
                Variations availabe:
                -'default'           -> Johnson-Champoux
                -'allard'            -> Johnson-Champoux-Allard
                -'lafarge'           -> Johnson-Champoux-Allard-Lafarge


        Returns:
        -------
            zc : int | float | complex
                Material Charactheristic Impedance
            kc : int | float | complex
                Material Wave Number
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

    elif var == 'allard':
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

    elif var == "lafarge":
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
