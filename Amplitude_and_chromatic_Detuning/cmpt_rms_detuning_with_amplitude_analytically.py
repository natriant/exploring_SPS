#Compute the rms detuning with amplitude from Eq. C6 and C7 of thesis

import numpy as np


def rms_amplitude_detuning_y_new(ey_geom, ex_geom, ayy, axy): # Eq.C7 my thesis
    DQ_rms_temp = 2*np.sqrt((ayy*ey_geom)**2+(axy*ex_geom)**2)
    return DQ_rms_temp



gamma_0 = 287.7 # 287.8 for 270 GeV
beta_0 = np.sqrt(1 - 1/gamma_0**2)

ex_norm_init = 2e-6 # m
ey_norm_init = 2e-6 # m


ex_geom, ey_geom = ex_norm_init/(beta_0*gamma_0), ey_norm_init/(beta_0*gamma_0)
print(ex_geom, ey_geom)

# klod =4
#ayy = 2782.5 # 1/m 
#axy = -1890.82

# klod = 1
ayy = 658.5
axy  = -769.43

Dqy_rms = rms_amplitude_detuning_y_new(ey_geom, ex_geom, ayy, axy)

print(Dqy_rms)
