import numpy as np
from matplotlib import pyplot as plt
import absorption_models as am  # Must be in the same folder as this file

# Setting up the variables to calculate the air properties
temp = 25  # Air temperature
hum = 30  # Air humidity
atm = 101320  # Atmospheric pressure

# Calculate the air properties
sound_speed, air_dens, c_imp_air, viscos, expans, Prandtl, \
    Cp = am.air_properties(temp, hum, atm)

# Set's up the frequency range that is desired
freq = np.arange(100, 10001, 1)

flow_resist = 35000  # Resistivity of the material
poros = 0.65  # Porosity of the material
tortu = 1  # Tortuosity of the material
visc = 500 * 10e-6  # Viscous characteristic length
term = 500 * 10e-6  # Thermal characteristic length
d = 0.05  # Material Thickness

# Johnson-Champoux formulation
zc, Kc = am.johnson_champoux(flow_resist, air_dens, poros, tortu, expans,
                             Prandtl, atm, visc, term, viscos, var='default')
absorption0 = am.johnson_champoux_absorption(zc, Kc, d, c_imp_air, poros)

# Johnson-Champoux-Allard formulation
zc, Kc = am.johnson_champoux(flow_resist, air_dens, poros, tortu, expans,
                             Prandtl, atm, visc, term, viscos, Cp, var='allard')
absorption1 = am.johnson_champoux_absorption(zc, Kc, d, c_imp_air, poros)

# Johnson-Champoux-Allard-Lafarge formulation
zc, Kc = am.johnson_champoux(flow_resist, air_dens, poros, tortu, expans,
                             Prandtl, atm, visc, term, viscos, Cp, var='lafarge')
absorption2 = am.johnson_champoux_absorption(zc, Kc, d, c_imp_air, poros)

# Putting all together
absorption = np.empty((absorption0.size, 3), dtype="complex64")
absorption[:, 0] = absorption0
absorption[:, 1] = absorption1
absorption[:, 2] = absorption2

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
