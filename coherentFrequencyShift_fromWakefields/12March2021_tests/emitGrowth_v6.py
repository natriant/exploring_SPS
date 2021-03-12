import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad
from scipy import constants as cst
import matplotlib
matplotlib.rcParams['font.size'] = 20

def hmm_gauss(freq,tau,m=0):
    return (2.0*np.pi*freq*tau/4)**(2*m) * np.exp(-(2.0*np.pi*freq*tau/4.)**2)

def hmm_gauss_omega(omega, tau, m=0):
    return (omega*tau/4)**(2*m) * np.exp(-(omega*tau/4.)**2)

if __name__ == '__main__':
    N = 3.5E10
    Qx, nu_b = 26.18, 0.18
    Qs = 0.0051

    circum = 2*np.pi*1.1E3
    frev = cst.c/circum
    omega_0 = 2*np.pi*frev

    gamma_t = 22.8
    energy = 270.0
    gamma = np.sqrt(1+(energy/0.938)**2)



    eta = 1/gamma_t**2 - 1/gamma**2
    tau = 1.85E-9
    sigmaz = tau*cst.c/4
    n_b, n, m_t = 1, 0, 1


    #### Impedance model from https://gitlab.cern.ch/IRIS/SPS_IW_model/-/tree/master/SPS_IW_model_python ###########################################
    impedanceData = np.genfromtxt('SPS_Complete_imp_model_2018_Q26.txt',skip_header=1,dtype=complex)
    freqZ = np.real(impedanceData[:, 0])*1E9  # frequencies in [GHz], so probably this needs to change in 1e9
    omegaZ = 2*np.pi*freqZ

    ReZ = np.real(impedanceData[:, 2])  # dipolar kick in y-plane
    ImZ = np.real(impedanceData[:, 2])
    #### Sacherer formula (e.g. 8 and 9 in Sec 2.5.7 in Handbook of Accelerator physics and engineering from A. Chao and M. Tigner) ################
    modeNumber = 0
    #nSideband = int(np.floor((1E10/frev)))

    nSideband = int(np.floor((1E10*2*np.pi / omega_0)))

    sidebands = np.arange(-nSideband, nSideband+0.5)

    chroma = 0.0
    chromaShift = chroma*omega_0/eta
    print(chromaShift)

    omega_mp = omega_0*(n_b*sidebands+n-m_t*Qx+modeNumber*Qs)
    #freqs = frev*(-0.18+sidebands+modeNumber*Qs)
    #hs = hmm_gauss(freqs-chromaShift, tau, m=modeNumber)
    hs = hmm_gauss_omega(omega_mp - chromaShift, tau, m=modeNumber)

    zeffs_dampingRate = np.interp(np.abs(omega_mp), omegaZ, ReZ)*np.sign(omega_mp)*hs
    zeffs_dampingRate /= np.sum(hs)
    zeff_dampingRate = np.sum(zeffs_dampingRate)

    zeffs_coherentDQ = np.interp(np.abs(omega_mp), omegaZ, ImZ) * hs
    zeffs_coherentDQ /= np.sum(hs)
    zeff_coherentDQ = np.sum(zeffs_coherentDQ)

    zeff = complex(zeff_dampingRate, zeff_coherentDQ)
    print(zeff)

    #print(f'zeff damping rate {zeff_dampingRate}')
    #print(f'zeff for coherent DQ {zeff_coherentDQ}')

    A = complex(0, - (cst.e ** 2 * N / (16.0 * np.pi * cst.m_p * gamma * Qx * frev * sigmaz * 2 * np.pi)))
    print(A)
    DQ = A * zeff
    print(DQ)




    '''
    #### Eq. 26 in https://aip.scitation.org/doi/abs/10.1063/1.47298 #################################################################################
    dGain = 2*dampingRate
    dmu = np.arange(1E-6,1E-5,5E-7)
    supps = np.zeros_like(dmu)
    for i in range(len(dmu)):
        f = lambda x : (4*np.pi**2*(1-dGain/2)**2*x**2)*np.exp(-x**2/(2.0*dmu[i]**2))/((4*np.pi**2*(1-dGain/2)*x**2+(dGain/2)**2)*np.sqrt(2*np.pi)*dmu[i])
        integral = quad(f,-10*dmu[i],10*dmu[i])
        supps[i] = integral[0]
    ##################################################################################################################################################

    fig = plt.figure(1)
    plt.plot(dmu*1E6,supps,'-b')
    plt.xlabel(r'R.m.s. tune spread [$10^{-6}$]')
    plt.ylabel(r'Emit. growth suppression factor')

    plt.show()
    '''