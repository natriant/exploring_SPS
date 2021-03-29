'''
Compute coherent tune shift from impedance using Eq.19 https://cds.cern.ch/record/322645/files/HEACC74_368-372.pdf from Sacherer's paper.
MKS (SI) units are used.
'''

import numpy as np
from scipy.constants import m_p, c, e
import matplotlib.pyplot as plt


def hmm_gaus(omega, sigma_z, l=0):
    return (omega*sigma_z/c)**(2*l)*np.exp(-(omega*sigma_z/c)**2)


#def hmm_gaus_tau(omega, tau):
#    return np.exp(-(omega*tau)**2)/np.sqrt(np.pi)


if __name__ == '__main__':
    l = 0 # azimuthial mode (headtail mode)
    circum = 2 * np.pi * 1.1E3  # [m[
    f_0 = c / circum  # revolution frequency in Hz
    omega_0 = 2 * np.pi * f_0  # angular revolution frequency


    # Load the impedance model

    #### Impedance model from https://gitlab.cern.ch/IRIS/SPS_IW_model/-/tree/master/SPS_IW_model_python ###########################################

    impedanceData = np.genfromtxt('SPS_Complete_imp_model_2018_Q26.txt', skip_header=1, dtype=complex)
    freqZ = np.real(impedanceData[:, 0]) * 1E9  # frequencies in [GHz], so probably this needs to change in 1e9
    ReZ = np.real(impedanceData[:, 2])  # dipolar kick in y-plane, always odd, -f(x) = f(-x)
    ImZ = np.imag(impedanceData[:, 2])  # dipolar kick in the y-plane, always even (f-x) = f(x)


    # A1. Plot vertical impedances (dipolar kick only)
    plt.plot(freqZ, ImZ, label=r'$\mathrm{Im(Z_y)}$')
    plt.plot(freqZ, ReZ, label=r'$\mathrm{Re(Z_y)}$')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()


    omegaZ = 2*np.pi*freqZ
    # A2. Plot vertical impedances (dipolar kick only)
    plt.plot(omegaZ, ImZ, label=r'$\mathrm{Im(Z_y)}$')
    plt.plot(omegaZ, ReZ, label=r'$\mathrm{Re(Z_y)}$')
    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()


    nSideband = 250000 #int(np.floor((1E10/(f_0))))
    print(nSideband)
    sidebands_p = np.arange(-nSideband, nSideband+0.5)


    # compute effective impedance
    #q_y = 0.18  # fractional part of the tune
    Q_y = 26.18  # betatron tune

    # bunch length
    #tau = 1.85E-9 # 4 sigma_t
    sigma_z = 0.155  # rms bunch length in [m] tau * c / 4
    tau = 4*sigma_z/c


    omegas = omega_0*(sidebands_p+Q_y) # the middle is not zero due to the shift of q_y
    #Qs =  0.0051
    #omegas = omega_0 * (sidebands_p + Q_y-l*Qs)


    Qp_y = 0.0


    gamma_t = 22.8  # for Q26

    # Compute lorentz facotr
    m0 = 0.9382720813e9 # proton rest mass, [eV/c^2]
    E_rest = m0
    E_0 = 270e9  # total energy of the reference particle [eV]
    gamma = E_0/E_rest

    eta = 1 / gamma_t ** 2 - 1 / gamma ** 2 # slip factor

    omega_xi = Qp_y * omega_0 / eta


    hs = hmm_gaus(omegas-omega_xi, sigma_z)
    #hs = hmm_gaus_tau(omegas-omega_xi, tau)


    plt.plot(omegas-omega_xi, hs)
    plt.show()
    plt.close()


    # You need to extrapolate Z eff in the omegas and also extend to negative frequencies
    #zeffs_coherentDQ = np.interp(np.abs(omega_mp), omegaZ, ImZ)

    # ImZ is always even f(x) = f(-x)
    ImZ_pos = ImZ
    ImZ_neg = ImZ


    # B. Plot vertical impedance also for negative frequencies
    plt.plot(omegaZ, ImZ_pos, label=r'$\mathrm{Im(Z_y)}$')
    plt.plot(-omegaZ, ImZ_neg, label=r'$\mathrm{-Im(Z_y)}$')

    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()
    plt.close()

    omegas_pos = list(filter(lambda x: x >= 0, omegas))
    omegas_neg = list(filter(lambda x: x < 0, omegas))

    ImZ_pos_interp = np.interp(omegas_pos, omegaZ, ImZ_pos)
    ImZ_neg_interp = np.interp(np.abs(omegas_neg), omegaZ, ImZ_neg)

    # C1. Plot ImZ(my_omegas)
    plt.plot(omegas_pos, ImZ_pos_interp, label=r'$\mathrm{Im(Z_y)}$')
    plt.plot(omegas_neg, ImZ_neg_interp, label=r'$\mathrm{-Im(Z_y)}$')

    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()
    plt.close()

    # C2. constract and plot total impedance.

    ImZ_interp_total = list(ImZ_neg_interp)+list(ImZ_pos_interp)

    plt.plot(omegas, ImZ_interp_total, label=r'$\mathrm{Im(Z_y)}$')
    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()
    plt.close()

    # compute effective impedance
    Zeff_nominator = np.sum(ImZ_interp_total*hs)
    Zeff_denominator = np.sum(hs)
    Zeff = Zeff_nominator/Zeff_denominator

    print(f'Zeff = {Zeff} [Ohm/m]')


    # Compute DQ
    intensity_list = np.linspace(0, 5e10, 5)
    Nb = 3.5e10 #intensity_list[4] #3e10  # protons per bunch
    I_0 = Nb*e*omega_0/(2*np.pi)  # bunch current
    beta = np.sqrt(1-(1/(gamma**2)))

    ### Eq.18 in https://cds.cern.ch/record/322645/files/HEACC74_368-372.pdf ###########
    Domega= -(e*beta*I_0*Zeff)/((1+l)*(2*26.18*gamma*m_p*4*sigma_z*omega_0))

    DQ_coh = Domega/omega_0
    print(f'DQ_coh = {DQ_coh} ')
    #DQ = -(beta*e*I_0*Zeff)/(4*sigma_z*np.sqrt(np.pi)*omega_0**2*gamma*26.18*m_p)
