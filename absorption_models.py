"""
Project description: a script containing all implemantations of the scientific
methods to found the absorption coefficient.

Project content:
- Delany Bazley (and with Miki variation and Allard-Champoux variation)
- Biot Allard
- Johnson-Champoux (and Allard variation)
- Rayleigh
"""

import numpy as np
from scipy import special as ss

formats = {"circle": 1,
           "square": 1.07,
           "equi-tri": 1.11,
           "retang": 0.81}


def air_properties(temp, hum, atm):
    KAPPA = 0.026  # W/(mK) air
    AIR_CONST = 287.031  # J/K/kg
    WATER_CONST = 461.521  # J/K/kg
    temp += 273.16  # K#
    pierce = 0.0658 * temp**3 - 53.7558 * temp**2 \
        + 14703.8127 * temp - 1345485.0465
    viscos = 7.72488e-8 * temp - 5.95238e-11 * temp**2 \
        + 2.71368e-14 * temp**3
    Cp = 4168.8 * (0.249679 - 7.55179e-5 * temp
                   + 1.69194e-7 * temp**2 - 6.46128e-11 * temp**3)
    # Volume Specific Heat Constant (J/kg/K) for 260 K < T < 600 K
    vol_Const = Cp - AIR_CONST
    Prandtl = viscos * Cp / KAPPA  # Prandtl number
    expans = Cp / vol_Const  # k
    airDens = atm / (AIR_CONST * temp) - (1 / AIR_CONST - 1 / WATER_CONST) * \
        hum / 100 * pierce / temp  # Air density
    soundSpd = (expans * atm / airDens)**0.5  # Speed of the sound
    characImpAir = soundSpd * airDens  # Characteristic Impedance of the Air
    return soundSpd, airDens, characImpAir, viscos, expans, Prandtl, Cp


def delany_bazley(freq, fluxRes, airDens, soundSpd, var='default'):
    if var == 'default':  # Original Delany Bazley
        R = 1 + 9.08 * ((1e3 * freq / fluxRes) ** -0.75)
        X = -11.9 * ((1e3 * freq / fluxRes) ** -0.73)
        Zc = soundSpd * airDens * (R + 1j * X)

        alpha = 1 + 10.8 * (1e3 * freq / fluxRes) ** -0.7
        beta = -10.3 * (1e3 * freq / fluxRes) ** -0.59
        kc = (2 * np.pi * freq / soundSpd) * (alpha + 1j * beta)

    elif var == 'miki':  # Miki variation
        R = 1 + 5.50 * ((1e3 * freq / fluxRes) ** -0.632)
        X = -8.43 * ((1e3 * freq / fluxRes) ** -0.632)
        Zc = soundSpd * airDens * (R + 1j * X)

        alpha = 1 + 7.81 * (1e3 * freq / fluxRes) ** -0.618
        beta = -11.41 * (1e3 * freq / fluxRes) ** -0.618
        kc = (2 * np.pi * freq / soundSpd) * (alpha + 1j * beta)

    elif var == 'allard-champoux':  # Allard and Champoux variation
        R = 1 + 0.0571 * (((airDens*freq) / fluxRes) ** -0.754)
        X = -0.0870 * (((airDens*freq) / fluxRes) ** -0.732)
        Zc = soundSpd * airDens * (R + 1j * X)

        alpha = 1 + 0.0978 * ((airDens*freq) / fluxRes) ** -0.7
        beta = -0.1890 * ((airDens*freq) / fluxRes) ** -0.595
        kc = (2 * np.pi * freq / soundSpd) * (alpha + 1j * beta)

    return Zc, kc


def delany_bazley_absorption(Zc, kc, d, cImpAir):
    Zs = -1j * (Zc / np.tan(kc * d))
    reflex = (Zs - cImpAir) / (Zs + cImpAir)
    absorption = 1 - np.abs(reflex) ** 2

    return absorption


def rayleigh(freq, fluxRes, airDens, soundSpd, poros):
    omega = 2 * np.pi * freq
    alpha = (1 - (1j * poros * fluxRes) / (airDens * omega)) ** 0.5
    kc = (omega/soundSpd) * alpha
    Zc = ((airDens * soundSpd)/poros) * alpha
    return Zc, kc


def rayleigh_absorption(Zc, kc, d, cImpAir):
    Zs = -1j * (Zc / np.tan(kc * d))
    reflex = (Zs - cImpAir) / (Zs + cImpAir)
    absorption = 1 - np.abs(reflex) ** 2
    return absorption


def shear_wave(omega, fluxRes, poros, tortus, shape, airDens):
    c1 = formats[shape]
    num = 8 * omega * airDens * tortus
    den = fluxRes * poros
    s = c1 * (num / den) ** 0.5

    return s


def biot_allard(freq, fluxRes, airDens, poros,
                tortus, expans, Prandtl, atm, shape):

    omega = 2 * np.pi * freq
    B = Prandtl ** 0.5
    s = shear_wave(omega, fluxRes, poros, tortus, shape, airDens)

    rhoEA = airDens * tortus
    rhoEB = (fluxRes * poros) / (1j * omega * airDens * tortus)
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

    Zc = (kE * rhoE) ** 0.5
    kc = omega * (rhoE / kE) ** 0.5

    return Zc, kc


def biot_allard_absorption(Zc, kc, d, cImpAir, poros):
    Zs = -1j * ((Zc/poros) / np.tan(kc * d))
    reflex = (Zs - cImpAir) / (Zs + cImpAir)
    absorption = 1 - np.abs(reflex) ** 2
    return absorption


def johnson_champoux(freq, fluxRes, airDens, poros, tort, expans, Prandtl, atm,
                     visc, term, neta, Cp=0, var='default'):
    KAPPA = 0.026  # W/(mK) air
    omega = 2 * np.pi * freq
    if var == 'default':
        # RhoE
        alpha = (1 - 1j * (4 * airDens * (tort**2) * neta * omega)
                 / ((fluxRes**2) * (poros**2) * (visc**2))) ** 0.5
        beta = 1 + (1j * (poros * fluxRes) / (airDens * omega * tort)) * alpha
        rhoE = airDens * tort * beta

        # kE
        epsilon = (1 - 1j * ((4 * airDens * Prandtl * (tort**2) * neta * omega)
                             / ((poros**2) * (fluxRes**2) * (term**2)))) ** 0.5
        gama = (airDens * omega * tort * Prandtl)
        zeta = (1 + (1j * ((poros * fluxRes) / gama)) * epsilon) ** -1
        eta = expans - (expans - 1) * zeta
        kE = (expans * atm) / eta

        # kc e zc
        kc = omega * ((rhoE / kE) ** 0.5)
        zc = (kE * rhoE) ** 0.5

        return zc, kc

    elif var == 'allard':
        # RhoE
        alpha = (1 + 1j * (4 * airDens * (tort**2) * neta * omega) /
                 ((fluxRes**2) * (poros**2) * (visc**2))) ** 0.5
        beta = 1 + (1j * (poros * fluxRes) / (airDens * omega * tort)) * alpha
        rhoE = ((airDens*tort)/poros) * beta

        # kE
        epsilon = (1 + 1j * ((airDens * (term**2) * Cp * omega) /
                             (16 *  KAPPA)) ** 0.5
        gama = (airDens * omega * (term**2) * Cp)
        zeta = (1 - (1j * ((8 * KAPPA) / gama)) * epsilon) ** -1
        eta = expans - (expans - 1) * zeta
        kE = ((expans * atm)/poros) / eta

        # kc e zc
        kc = omega * ((rhoE / kE) ** 0.5)
        zc = (kE * rhoE) ** 0.5

        return zc, kc


def johnson_champoux_absorption(zc, kc, cImpAir, poros, d):
    zs = 1j * (zc / poros) * (1 / np.tan(kc * d))
    reflex = (zs - cImpAir) / (zs + cImpAir)
    absorption = 1 - np.abs(reflex) ** 2
    return absorption
