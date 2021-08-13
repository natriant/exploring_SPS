import sys
sys.path.append('../../utils/')
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from scipy.constants import c
from math import * 

from coordinatesConversions import *
from cmpt_TuneSpreads import *

### Plotting parameters #####
params = {'legend.fontsize': 30,
          'figure.figsize': (9.5, 8.5),
          'axes.labelsize': 30,
          'axes.titlesize': 30,
          'xtick.labelsize': 30,
          'ytick.labelsize': 30,
          'image.cmap': 'jet',
          'lines.linewidth': 5,
          'lines.markersize': 12,
          'font.family': 'sans-serif'}

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rcParams.update(params)

#### my functions #######3
def cmpt_Brho(pc):
    # pc in GeV
    return pc/(c*1e-9) # from madx


def cmpt_octupole_coefficient(k3, Brho):
    return k3*Brho

def cmpt_current_klof(O3): # O3 the octupole coefficeintt
    return 500*O3/11779.66

def cmpt_current_klod(O3): # O3 the octupole coefficeintt
    return 500*O3/44490.4


#### Study parameters #####
pc = 270 # GeV
Brho = cmpt_Brho(pc)

study_list = ['matching_results_constraintaxy0_nobb5b7_noChroma.pkl', 'matching_results_NOconstraintaxy_nobb5b7_noChroma.pkl', 'matching_results_constraintaxy0_nobb5b7_noChroma_negativeayy.pkl', 'matching_results_Noconstraintaxy_nobb5b7_noChroma_negativeayy.pkl']

fig, ax = plt.subplots(1,1)
ax2 = ax.twiny()

### for legend ###
ax.scatter(0, 0, color = 'g', label=r'$\mathrm{I_{klof} \ [A] \leq 100} $')
ax.scatter(0, 0, color = 'y', label=r'$\mathrm{100 < I_{klof} \ [A] \leq 200} $')
ax.scatter(0, 0, color = 'r', label=r'$\mathrm{200  < I_{klof} \ [A]} $')

# Load data from octupole matching
for study in study_list:
    df = pd.read_pickle(study)
    my_x, my_y = df.ayy, df.axy

    ## Compute DQrms from amplitude detuning ###
    # A. Load a typical SPS bunch used in the MD
    path_to_bunch = '/home/natalia/PhD_projects/cc_rf_noise_pyheadtail/post_process_analysis/studies_with_sets_of_kicks/'
    bunch = pickle.load(open(path_to_bunch+'bunch_IPAC', 'rb'))

    # Optics at CC2
    beta_y = 73.81671646
    beta_x = 30.31164764
    alpha_y = 0
    alpha_x = 0

    # Coordinates
    x, px = bunch.x, bunch.xp
    y, py =  bunch.y, bunch.yp

    # Normalised coordinates 
    x_n, px_n = cmpt_normalised_coordinates(x, px, beta_x, alpha_x)
    y_n, py_n = cmpt_normalised_coordinates(y, py, beta_y, alpha_y)


    # Compute actions
    Jx_init = cmpt_actions(x_n, px_n)
    Jy_init = cmpt_actions(y_n, py_n)

    rms_Jx_init = np.std(Jx_init)
    rms_Jy_init = np.std(Jy_init)

    ############
    klof, klod = df.klof, df.klod

    O3_klof = cmpt_octupole_coefficient(klof, Brho)
    O3_klod = cmpt_octupole_coefficient(klod, Brho)

    I_klof = list(cmpt_current_klof(O3_klof))
    I_klod = list(cmpt_current_klod(O3_klod))

    colors = []
    Dqy_rms = []

    for i, my_i_klof in enumerate(I_klof):
        Dqy_rms.append(rms_amplitude_detuning_y(Jy_init, Jx_init, my_y[i], my_x[i]))
        if np.abs(my_i_klof) <=100:
            colors.append('g')
        elif np.abs(my_i_klof) > 100 and np.abs(my_i_klof) <=200:
            colors.append('y')
        else:
            colors.append('r')

    for i in range(len(I_klof)):
        ax.scatter(my_x[i], my_y[i], color=colors[i])

    if study == study_list[0]:
        my_dqy = Dqy_rms
        my_jy_rms, my_jx_rms = rms_Jy_init, rms_Jx_init
        print(Dqy_rms)
        ax2.plot(np.array(Dqy_rms)*1e4, my_y, linestyle='')
        #ax2.plot(-np.array(Dqy_rms)*1e4, np.zeros(len(Dqy_rms)), linestyle='')
    

ax2.set_xticks(np.arange(-3, 4,1))
ax2.set_xticklabels(np.abs((np.arange(-3, 4,1))))

#ax2.set_xlim(-np.max(my_dqy)-21000*2*my_jy_rms*1e4, np.max(my_dqy)+21000*2*my_jy_rms*1e4 )
#ax.set_xlim(-21000, 21000)

ax.set_xlabel(r'$\mathrm{\alpha_{yy}}$'+' [1/m]')
ax.set_ylabel(r'$\mathrm{\alpha_{xy}}$'+' [1/m]')
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax.grid(ls='--')
ax.legend(frameon=False)
plt.tight_layout()
plt.savefig('possible_octuople_settings_klof.png')

plt.show()




