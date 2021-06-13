import numpy as np
from matplotlib import pyplot as plt
import pyabsorp as ab  # Must be in the same folder as this file

# Setting up the variables to calculate the air properties
temp = 25  # Air temperature
hum = 30  # Air humidity
atm = 101320  # Atmospheric pressure

# Calculate the air properties
sound_speed, air_dens, z_air, viscos, gama, Prandtl, \
    Cp = ab.air_properties(temp, hum, atm)

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
zc, Kc = ab.johnson_champoux(flow_resist, air_dens, poros, tortu, gama,
                             Prandtl, atm, visc, term, viscos, var='default')
absorption0 = ab.absorption_coefficient(zc, Kc, d, z_air)

# Johnson-Champoux-Allard formulation
zc, Kc = ab.johnson_champoux(flow_resist, air_dens, poros, tortu, gama,
                             Prandtl, atm, visc, term, viscos, Cp, var='allard')
absorption1 = ab.absorption_coefficient(zc, Kc, d, z_air)

# Johnson-Champoux-Allard-Lafarge formulation
zc, Kc = ab.johnson_champoux(flow_resist, air_dens, poros, tortu, gama,
                             Prandtl, atm, visc, term, viscos, therm_perm, 
                             Cp, var='lafarge')
absorption2 = ab.absorption_coefficient(zc, Kc, d, z_air)

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
