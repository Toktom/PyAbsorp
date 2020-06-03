"""
Author: Michael Markus Ackermann
================================
Here you will find everything related to the johnson-champoux model.
"""


import numpy as np


def johnson_champoux(flow_resis, air_dens, poros, tortu, expans,
                     prandtl, atm, visc, term, neta, Cp=0,
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
            expans: int | float
                Ratio of specific heat
            prandtl: int | float
                Prandtl's number
            atm: int
                Atmospheric pressure
            visc: int | float
                Viscous characteristic length
            term: int | float
                Thermal characteristic length
            neta: int | float
                Dynamic vicosity of air
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
    KAPPA = 0.026  # W/(mK) air
    STATIC_THERM_PERM = 23e-10  # Static thermal permeability
    omega = 2 * np.pi * freq

    if var == 'default':

        gama = 4 * omega * air_dens * neta * (tortu**2)

        # rho_ef
        alpha = (1 - 1j * gama
                 / ((flow_resis**2) * (poros**2) * (visc**2))) ** 0.5
        delta = air_dens * omega * tortu
        beta = 1 - (1j * (poros * flow_resis) / delta) * alpha
        rho_ef = air_dens * tortu * beta

        # k_ef
        epsilon = (1 + 1j * ((gama * prandtl)
                   / ((poros**2) * (flow_resis**2) * (term**2)))) ** 0.5
        psi = air_dens * omega * tortu * prandtl
        zeta = (1 - (1j * ((poros * flow_resis) / psi) * epsilon)) ** -1
        eta = expans - (expans - 1) * zeta
        k_ef = (expans * atm) / eta

        # Changing from efficient to equivalent
        rho_eq = rho_ef / poros
        k_eq = k_ef / poros

        # Charactheristic Impedance (zc) and the Wave Number (kc)
        kc = omega * ((rho_eq / k_eq) ** 0.5)
        zc = (k_eq * rho_eq) ** 0.5

        return zc, kc

    elif var == 'allard':
        # rho_ef
        alpha = (1 - 1j * (4 * air_dens * (tortu**2) * neta * omega)
                 / ((flow_resis**2) * (poros**2) * (visc**2))) ** 0.5
        beta = 1 - (1j * (poros * flow_resis) /
                    (air_dens * omega * tortu)) * alpha
        rho_ef = ((air_dens*tortu) / poros) * beta

        # k_ef
        epsilon = (1 + 1j * ((air_dens * (term**2) * Cp * omega)
                   / (16 * KAPPA))) ** 0.5
        gama = (air_dens * omega * (term**2) * Cp)
        zeta = (1 - (1j * ((8 * KAPPA) / gama)) * epsilon) ** -1
        eta = expans - (expans - 1) * zeta
        k_ef = ((expans * atm)/poros) / eta

        # Changing from efficient to equivalent
        rho_eq = rho_ef / poros
        k_eq = k_ef / poros

        # Charactheristic Impedance (zc) and the Wave Number (kc)
        kc = omega * ((rho_eq / k_eq) ** 0.5)
        zc = (k_eq * rho_eq) ** 0.5

        return zc, kc

    elif var == "lafarge":
        # rho_ef
        alpha = (1 - 1j * (4 * air_dens * (tortu**2) * neta * omega)
                 / ((flow_resis**2) * (poros**2) * (visc**2))) ** 0.5
        beta = 1 - (1j * (poros * flow_resis) /
                    (air_dens * omega * tortu)) * alpha
        rho_ef = ((air_dens*tortu) / poros) * beta

        # k_ef
        psi = 4 * (STATIC_THERM_PERM**2) * Cp * air_dens * omega
        epsilon = (1 - 1j * ((psi) / (KAPPA * (term**2) * (poros**2)))) ** 0.5
        gama = (STATIC_THERM_PERM * air_dens * Cp * omega)
        zeta = (1 + (1j * ((poros * KAPPA) / gama)) * epsilon) ** -1
        eta = expans - (expans - 1) * zeta
        k_ef = ((expans * atm)/poros) / eta

        # Changing from efficient to equivalent
        rho_eq = rho_ef / poros
        k_eq = k_ef / poros

        # Charactheristic Impedance (zc) and the Wave Number (kc)
        kc = omega * ((rho_eq / k_eq) ** 0.5)
        zc = (k_eq * rho_eq) ** 0.5

        return zc, kc
