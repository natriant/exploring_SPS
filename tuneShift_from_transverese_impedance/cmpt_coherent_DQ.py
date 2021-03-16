import numpy as np
from scipy.constants import m_p, c, e
import matplotlib.pyplot as plt

Nb = 3e10 # protons per bunch
q_y = 0.18 # fractional part of the tune

circum = 2 * np.pi * 1.1E3  # [m[
f_0 = c / circum  # revolution frequency in Hz
omega_0 = 2 * np.pi * f_0  # angular revolution frequency

I_0 = Nb*e*omega_0/(2*np.pi)  # bunch current

#tau = 1.85E-9 # 4 sigma_t
sigma_z = 0.155  # rms bunch length in [m] tau * c / 4


hmm_gauss = lambda omega: e**(-(omega*sigma_z)**2/(c**2))  #call by: hmm_gaus(omega_values)

my_omegas = lambda p: omega_0*(p+q_y)  # the discrete frequencies of ocsillations

# Load the impedance model

#### Impedance model from https://gitlab.cern.ch/IRIS/SPS_IW_model/-/tree/master/SPS_IW_model_python ###########################################

impedanceData = np.genfromtxt('SPS_Complete_imp_model_2018_Q26.txt', skip_header=1, dtype=complex)
freqZ = np.real(impedanceData[:, 0]) * 1E9  # frequencies in [GHz], so probably this needs to change in 1e9
ReZ = np.real(impedanceData[:, 2])  # dipolar kick in y-plane
ImZ = np.imag(impedanceData[:, 2])

# Plot vertical impedances (dipolar kick only)
plt.plot(freqZ, ImZ, label=r'$\mathrm{Im(Z_y)}$')
plt.plot(freqZ, ReZ, label=r'$\mathrm{Re(Z_y)}$')
plt.xlabel('Frequency [Hz]')
plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
plt.legend()
plt.grid()
plt.ylim(-0.5e7, 0.5e8)
plt.show()


