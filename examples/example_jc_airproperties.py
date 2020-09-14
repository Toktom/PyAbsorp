#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to demonstrate how to use the AirProperties class within
the same example of the Johnson-Champoux absorption model.

Created on Thu Jul 30 00:47:12 2020

author: João Vitor Gutkoski Paes
"""

import numpy as np
from matplotlib import pyplot as plt
import pyabsorp as ab  # Must be in the same folder as this file

# Using new AirProperties class
air = ab.AirProperties(temp=25, hum=30, atm=101320)

# Set's up the frequency range that is desired
freq = np.arange(100, 10001, 1)

flow_resist = 35000  # Resistivity of the material
poros = 0.65  # Porosity of the material
tortu = 1  # Tortuosity of the material
visc = 750 * 10e-6  # Viscous characteristic length
term = 500 * 10e-6  # Thermal characteristic length
d = 0.05  # Material Thickness
therm_perm = poros * (term ** 2) / 8  # cylindrical pore

# Johnson-Champoux formulation
zc, Kc = ab.johnson_champoux(flow_resist, air.density, poros, tortu,
                             air.specHeatRatio, air.prandtl, air.atmPressure,
                             visc, term, air.viscosity, var='default')
absorption0 = ab.absorption_coefficient(zc, Kc, d, air.impedance)

# Johnson-Champoux-Allard formulation
zc, Kc = ab.johnson_champoux(flow_resist, air.density, poros, tortu,
                             air.specHeatRatio, air.prandtl, air.atmPressure,
                             visc, term, air.viscosity, air.specHeatCP, var='allard')
absorption1 = ab.absorption_coefficient(zc, Kc, d, air.impedance)

# Johnson-Champoux-Allard-Lafarge formulation
zc, Kc = ab.johnson_champoux(flow_resist, air.density, poros, tortu,
                             air.specHeatRatio, air.prandtl, air.atmPressure,
                             visc, term, air.viscosity, therm_perm, air.specHeatCP,
                             var='lafarge')
absorption2 = ab.absorption_coefficient(zc, Kc, d, air.impedance)

# Putting all together
absorption = np.hstack((absorption0[:, None], absorption1[:, None], absorption2[:, None]))

legends = ['\u03B1 Johnson-Chmapoux',
           '\u03B1 Johnson-Chmapoux-Allard',
           '\u03B1 Johnson-Chmapoux-Allard-Lafarge']

plt.title("Sound Absorption Coefficient Chart")
for idc, lgnd in enumerate(legends):
    plt.semilogx(freq, 100 * absorption[:, idc], label=lgnd)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Absorption Coefficient [%]')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.show()

