
"""
Project description: a script containing all implemantations of the scientific
methods to found the absortion coefficient.

Project content:
- Delany Bazley
- Biot Allard
- Johnson-Champoux
"""

import numpy as np
from scipy import special as ss

formats = {"circle": 1,
           "square": 1.07,
           "equi-tri": 1.11,
           "retang": 0.81}


def air_properties(temp, hum, atm):
    thermCond = 0.026  # W/(mK)
    temp += 273.16  # K
    airConst = 287.031  # J/K/kg
    waterConst = 461.521  # J/K/kg
    Pierce = 0.0658 * temp**3 - 53.7558 * temp**2 \
        + 14703.8127 * temp - 1345485.0465
    viscos = 7.72488e-8 * temp - 5.95238e-11 * temp**2 \
        + 2.71368e-14 * temp**3
    pressConst = 4168.8 * (0.249679 - 7.55179e-5 * temp
                           + 1.69194e-7 * temp**2 - 6.46128e-11 * temp**3)
    # Volume Specific Heat Constant (J/kg/K) for 260 K < T < 600 K
    volConst = pressConst - airConst
    Prandtl = viscos * pressConst / thermCond  # Prandtl number
    expans = pressConst / volConst  # k
    airDens = atm / (airConst * temp) - (1 / airConst - 1 / waterConst) * \
        hum / 100 * Pierce / temp  # Air density
    soundSpd = (expans * atm / airDens)**0.5  # Speed of the sound
    characImpAir = soundSpd * airDens  # Characteristic Impedance of the Air
    return soundSpd, airDens, characImpAir, viscos, expans, Prandtl, pressConst


def delany_bazley(freq, fluxRes, soundSpd, airDens):
    # Characteristic Impedance of the Material (Zc)
    # Zc = Zo ( R + jX )
    R = 1 + 9.08 * ((1e3 * freq / fluxRes) ** -0.75)
    X = -11.9 * ((1e3 * freq / fluxRes) ** -0.73)
    Zc = soundSpd * airDens * (R + 1j * X)

    # Complex wave number (kc)
    # kc = k0 ( Alpha + jBeta )
    Alpha = 1 + 10.8 * (1e3 * freq / fluxRes) ** -0.7
    Beta = -10.3 * (1e3 * freq / fluxRes) ** -0.59
    kc = (2 * np.pi * freq / soundSpd) * (Alpha + 1j * Beta)

    return Zc, kc


def delany_bazley_absortion(Zc, kc, d, cImpAir):
    Zs = -1j * (Zc / np.tan(kc * d))
    reflex = (Zs - cImpAir) / (Zs + cImpAir)
    absortion = 1 - np.abs(reflex) ** 2
    return absortion


def shear_wave(omega, fluxRes, poros, tortus, shape, airDensAr):
    c1 = formats[shape]
    num = 8 * omega * airDensAr * tortus
    den = fluxRes * poros
    s = c1 * (num / den)**0.5
    return s


def biot_allard(freq, fluxRes, poros, tortus, shape, airDens, expans, Prandtl, atm):

    omega = 2 * np.pi * freq
    B = Prandtl**0.5
    s = shear_wave(omega, fluxRes, poros, tortus, shape, airDens)

    airDensA = airDens * tortus
    airDensB = (fluxRes * poros) / (1j * omega * airDens * tortus)
    airDensC = (s * (-1j) ** 0.5) / 4
    airDensD = 2 / (s * (-1j)**0.5)
    airDensE = ss.jv(1, s * (-1j)**0.5) / ss.jv(0, s * (-1j)**0.5)

    airDensity = airDensA * (1 - airDensB * ((airDensC * airDensE)
                                            / (1 - airDensD * airDensE)))
    compressB = airDensB / B
    compressC = airDensC * B
    compressD = airDensD / B
    compressE = ss.jv(1, s * B * (-1j)**0.5) / ss.jv(0, s * B * (-1j)**0.5)

    compress = (expans * atm) \
        / (expans - (expans - 1) / (1 - compressB * ((compressC * compressE)
                                             / (1 - compressD * compressE))))

    Zc = (compress * airDensity)**0.5
    kc = omega * (airDensity / compress)**0.5
    return Zc, kc


def biot_allard_absortion(Zc, kc, d, cImpAir):
    Zs = -1j * (Zc / np.tan(kc * d))
    reflex = (Zs - cImpAir) / (Zs + cImpAir)
    absortion = 1 - np.abs(reflex) ** 2
    return absortion


def johnson_champoux(fluxRes, tort, poros, visc, term, airDens, neta, Prandtl, expans, atm, freq):
    omega = 2 * np.pi * freq  # Angular Frequency

    # RhoE
    alpha = (1 - 1j * (4 * airDens * (tort**2) * neta * omega) /
             ((fluxRes**2) * (poros**2) * (visc**2)))**(1 / 2)
    beta = 1 + (1j * (poros * fluxRes) / (airDens * omega * tort)) * alpha
    rhoE = airDens * tort * beta

    # kE
    epsilon = (1 - 1j * ((4 * airDens * Prandtl * (tort**2) * neta * omega) /
                             ((poros**2) * (fluxRes**2) * (term**2))))**(1 / 2)
    zeta = (1 + (1j * ((poros * fluxRes) /
                         (airDens * omega * tort * Prandtl))) * epsilon) ** -1
    eta = expans - (expans - 1) * zeta
    kE = (expans * atm) / eta

    # kc e zc
    kc = omega * ((rhoE / kE)**(1 / 2))
    zc = (kE * rhoE)**(1 / 2)

    return zc, kc


def johnson_champoux_absortion(Zc, kc, cImpAir, poros, d):
    Zs = 1j * (Zc / poros) * (1 / np.tan(kc * d))
    reflex = (Zs - cImpAir) / (Zs + cImpAir)
    absortion = 1 - np.abs(reflex) ** 2
    return absortion
