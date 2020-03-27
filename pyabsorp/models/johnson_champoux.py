"""

=================
@Author: Michael Markus Ackermann
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

        gama = 4 * air_dens * (tortu**2) * neta * omega

        # RhoE
        alpha = (1 - 1j * gama
                 / ((flow_resis**2) * (poros**2) * (visc**2))) ** 0.5
        delta = air_dens * omega * tortu
        beta = 1 + (1j * (poros * flow_resis) / delta) * alpha
        rhoE = air_dens * tortu * beta

        # kE
        epsilon = (1 - 1j * ((gama * prandtl)
                   / ((poros**2) * (flow_resis**2) * (term**2)))) ** 0.5
        psi = air_dens * omega * tortu * prandtl
        zeta = (1 + (1j * ((poros * flow_resis) / psi) * epsilon)) ** -1
        eta = expans - (expans - 1) * zeta
        kE = (expans * atm) / eta

        # kc e zc
        kc = omega * ((rhoE / kE) ** 0.5)
        zc = (kE * rhoE) ** 0.5

        return zc, kc

    elif var == 'allard':
        # RhoE
        alpha = (1 + 1j * (4 * air_dens * (tortu**2) * neta * omega)
                 / ((flow_resis**2) * (poros**2) * (visc**2))) ** 0.5
        beta = 1 + (1j * (poros * flow_resis) /
                    (air_dens * omega * tortu)) * alpha
        rhoE = ((air_dens*tortu) / poros) * beta

        # kE
        epsilon = (1 + 1j * ((air_dens * (term**2) * Cp * omega)
                   / (16 * KAPPA))) ** 0.5
        gama = (air_dens * omega * (term**2) * Cp)
        zeta = (1 - (1j * ((8 * KAPPA) / gama)) * epsilon) ** -1
        eta = expans - (expans - 1) * zeta
        kE = ((expans * atm)/poros) / eta

        # kc e zc
        kc = omega * ((rhoE / kE) ** 0.5)
        zc = (kE * rhoE) ** 0.5

        return zc, kc

    elif var == "lafarge":
        # RhoE
        alpha = (1 + 1j * (4 * air_dens * (tortu**2) * neta * omega)
                 / ((flow_resis**2) * (poros**2) * (visc**2))) ** 0.5
        beta = 1 + (1j * (poros * flow_resis) /
                    (air_dens * omega * tortu)) * alpha
        rhoE = ((air_dens*tortu) / poros) * beta

        # kE
        psi = 4 * (STATIC_THERM_PERM**2) * Cp * air_dens * omega
        epsilon = (1 + 1j * ((psi) / (KAPPA * (term**2) * (poros**2)))) ** 0.5
        gama = (STATIC_THERM_PERM * air_dens * Cp * omega)
        zeta = (1 - (1j * ((poros * KAPPA) / gama)) * epsilon) ** -1
        eta = expans - (expans - 1) * zeta
        kE = ((expans * atm)/poros) / eta

        # kc e zc
        kc = omega * ((rhoE / kE) ** 0.5)
        zc = (kE * rhoE) ** 0.5

        return zc, kc


def johnson_champoux_absorption(zc, kc, d, c_imp_air, poros):
    """
    Returns the Sound Absorption Coefficient for the Johnson-Champoux Model.
    NOTE: Only use it with the johnson_champoux function.

        Parameters:
        -----------
            zc : int | float | complex
                Material Charactheristic Impedance
            kc : int | float | complex
                Material Wave Number
            d : float
                Material Thickness
            c_imp_air : int | float
                Air Characteristic Impedance
            poros: float
                Porosity of the Material

        Returns:
        --------
            absorption : int | ndarray
                Sound Absorption Coefficient of the Material
    """
    zs = 1j * (zc / poros) * (1 / np.tan(kc * d))
    reflex = (zs - c_imp_air) / (zs + c_imp_air)
    absorption = 1 - np.abs(reflex) ** 2
    return absorption
