import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import c

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
### for legend ###
ax.scatter(0, 0, color = 'g', label=r'$\mathrm{I_{klod} \ and \  I_{klof} \ [A] \leq 100} $')
ax.scatter(0, 0, color = 'y', label=r'$\mathrm{100 < I_{klod} \ or \  I_{klof} \ [A] \leq 200} $')
ax.scatter(0, 0, color = 'r', label=r'$\mathrm{200  < I_{klod} \ or \ I_{klof} \ [A]} $')

# Load data from octupole matching
for study in study_list:
    df = pd.read_pickle(study)
    my_x, my_y = df.ayy, df.axy
    klof, klod = df.klof, df.klod

    O3_klof = cmpt_octupole_coefficient(klof, Brho)
    O3_klod = cmpt_octupole_coefficient(klod, Brho)

    I_klof = list(cmpt_current_klof(O3_klof))
    I_klod = list(cmpt_current_klod(O3_klod))

    colors = []
    for i, my_i_klod in enumerate(I_klod):
        my_i_klof = I_klof[i]
        if np.abs(my_i_klod) > 200 or np.abs(my_i_klof) > 200 :
            colors.append('r')
        elif (np.abs(my_i_klod) > 100 and np.abs(my_i_klod) <=200) or (np.abs(my_i_klof) > 100 and np.abs(my_i_klof) <=200):
            colors.append('y')
        else:
            colors.append('g')

    for i in range(len(I_klod)):
        ax.scatter(my_x[i], my_y[i], color=colors[i])


ax.set_xlabel(r'$\mathrm{\alpha_{yy}}$'+' [1/m]')
ax.set_ylabel(r'$\mathrm{\alpha_{xy}}$'+' [1/m]')
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax.grid(ls='--')
ax.legend(frameon=False, loc=3)
plt.tight_layout()
plt.savefig('possible_octuople_settings_klod_klof.png')

plt.show()




