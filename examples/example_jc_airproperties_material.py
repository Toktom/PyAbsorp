#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script to estimate a material absorption coeficient.
Using AirProperties and Material classes within same Johnson-Champoux example code.

Created on Thu Jul 30 00:55:08 2020

author: João Vitor Gutkoski Paes
"""

import numpy as np
from matplotlib import pyplot as plt
import pyabsorp as ab  # Must be in the same folder as this file


# Definition of the problem
# freq = np.arange(100, 10001, 1)          # Frequencies

# Using new AirProperties class
air = ab.Air(temperature=25, humidity=30, atmospheric_pressure=101320)

# Using new Material class
# Most of these values should be acquired from laboratory measurements.
# Bibliography values are also good for code testing and validation.
mat = ab.Material(air=air, thickness=0.05, flow_resistivity=35000, porosity=0.65,
                  tortuosity=1., viscosity_length=750e-5, thermal_length=500e-5,
                  thermal_permeability=((0.65/8)*(500e-5)**2))


variations = ['default', 'allard', 'lafarge']

legends = ['\u03B1 Johnson-Chmapoux',
           '\u03B1 Johnson-Chmapoux-Allard',
           '\u03B1 Johnson-Chmapoux-Allard-Lafarge']


plt.title("Sound Absorption Coefficient Chart")

for label, var in zip(legends, variations):
    # Calculate the absorption for a given variation
    mat.estimate_absorption(mat.frequencies, method='jc', var=var)

    # Draw the lines
    plt.semilogx(mat.frequencies, 100*mat.absorption, label=label)

plt.xlabel('Frequency [Hz]')
plt.ylabel('Absorption Coefficient [%]')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.show()
