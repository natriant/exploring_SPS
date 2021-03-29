''''
Compute damping/growth time from impedance using Eq.6.207  https://www.slac.stanford.edu/~achao/WileyBook/WileyChapter6.pdf of Chao book.
CGS units are used.
'''

import os
import pandas as pd
import numpy as np
from scipy.constants import m_p, c, e
import matplotlib.pyplot as plt
from scipy.integrate import quad



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


    nSideband = int(np.floor((1E10/(f_0))))
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


    Qp_y = 1.0


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

    # ReZ is always odd -f(x) = f(-x)
    ReZ_pos = ReZ
    ReZ_neg = -ReZ


    # B. Plot vertical impedance also for negative frequencies
    plt.plot(omegaZ, ReZ_pos, label=r'$\mathrm{Re(Z_y)}$')
    plt.plot(-omegaZ, ReZ_neg, label=r'$\mathrm{-Re(Z_y)}$')

    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()
    plt.close()


    omegas_pos = list(filter(lambda x: x >=0, omegas))
    omegas_neg = list(filter(lambda x: x < 0, omegas))

    ReZ_pos_interp = np.interp(omegas_pos, omegaZ, ReZ_pos)
    ReZ_neg_interp = np.interp(np.abs(omegas_neg), omegaZ, ReZ_neg)

    # C1. Plot ImZ(my_omegas)
    plt.plot(omegas_pos, ReZ_pos_interp, label=r'$\mathrm{Re(Z_y)}$')
    plt.plot(omegas_neg, ReZ_neg_interp, label=r'$\mathrm{-Re(Z_y)}$')

    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()
    plt.close()

    # C2. constract and plot total impedance.

    ReZ_interp_total = list(ReZ_neg_interp)+list(ReZ_pos_interp)

    plt.plot(omegas, ReZ_interp_total, label=r'$\mathrm{Im(Z_y)}$')
    plt.xlabel('Angular frequency [rad]')
    plt.ylabel(r'$\mathrm{Z_y \ [\Omega /m]}$')
    plt.legend()
    plt.grid(ls='--')
    plt.ylim(-0.5e7, 0.5e8)
    plt.show()
    plt.close()

    # compute effective impedance
    Zeff_nominator = np.sum(ReZ_interp_total*hs)
    Zeff_denominator = np.sum(hs)
    Zeff = Zeff_nominator/Zeff_denominator

    print(Zeff)


    # Compute DQ

    Nb = 3.5e10  # protons per bunch
    I_0 = Nb*e*omega_0/(2*np.pi)  # bunch current
    beta = np.sqrt(1-(1/(gamma**2)))

    ### Eq.18 in https://cds.cern.ch/record/322645/files/HEACC74_368-372.pdf ###########
    Domega= (e*beta*I_0*Zeff)/((1+l)*(2*26.18*gamma*m_p*4*sigma_z*omega_0))

    growth = -Domega/omega_0
    print(f'Growth rate 1/τ = {growth} ')
    #DQ = -(beta*e*I_0*Zeff)/(4*sigma_z*np.sqrt(np.pi)*omega_0**2*gamma*26.18*m_p)

    dGain = 2 * np.abs(growth)
    dmu = np.arange(1E-6, 3E-4, 1E-6)  # tune spread (striclty from amplitude detuning)
    supps = np.zeros_like(dmu)

    #### Eq. 26 in https://aip.scitation.org/doi/abs/10.1063/1.47298 ####################################

    #for i in range(len(dmu)):
    #    f = lambda x: (4 * np.pi ** 2 * (1 - dGain / 2) ** 2 * x ** 2) * np.exp(-x ** 2 / (2.0 * dmu[i] ** 2)) / (
    #                (4 * np.pi ** 2 * (1 - dGain / 2) * x ** 2 + (dGain / 2) ** 2) * np.sqrt(2 * np.pi) * dmu[i])
    #    integral = quad(f, -10 * dmu[i], 10 * dmu[i])
    #    supps[i] = integral[0]


    for i in range(len(dmu)):
        supps[i] = (4*np.pi**2*(1-dGain/2)**2*dmu[i]**2)/(4*np.pi**2*(1-dGain/2)*dmu[i]**2+(dGain/2)**2)
    ##################################################################################################################################################

    fig = plt.figure(1)
    #plt.hlines(1, 1, 2)
    plt.plot(dmu * 1E4, supps, '-b')
    plt.xlabel(r'R.m.s. tune spread [$10^{-4}$]')
    plt.ylabel(r'Emit. growth suppression factor')
    plt.tight_layout()
    plt.show()

    save_supps = False
    if save_supps:
        data = {'tune spread value': list(dmu), 'suppression factor':list(supps)}
        df = pd.DataFrame(data)
        print(df)
        filename = f'suppressionFactor_Intensity{Nb/1e10}e10_sigma_z{sigma_z}_Qpy{Qp_y}_dmuScan.pkl'
        if not os.path.exists(filename):
            df.to_pickle(filename)
