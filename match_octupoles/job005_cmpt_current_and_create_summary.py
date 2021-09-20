'''
Approach 1 and 2: https://www.overleaf.com/2266482692wykkztbjdzhm
'''
from cpymad.madx import Madx
import pandas as pd
import numpy as np
from scipy.constants import c as clight

#### my functions #######
def cmpt_Brho(pc):
    # pc in GeV
    return pc/(clight*1e-9) # from madx

def cmpt_octupole_coefficient(k3, Brho): # [T/m^3]
    return k3*Brho

def cmpt_current_klof(O3): # O3 the octupole coefficient [A]
    # return 500*O3/11779.66 # Approach 1
    return 0.037*O3 # Approach 2

def cmpt_current_klod(O3): # O3 the octupole coefficient [A]
    #return 500*O3/44490.4 # Approach 1
    return 0.01*O3 # Approach 2

#### Study parameters #####
pc = 270 # GeV
Brho = cmpt_Brho(pc)

path2octupoles_matching_data  = './'
study_list = ['matching_results_QpxQpy1_b3b5b7_270GeV_negative_ayy_axyNoConstraint_lod.pkl', 'matching_results_QpxQpy1_b3b5b7_270GeV_positive_ayy_axyNoConstraint_lod.pkl']

axy_list, ayy_list = [], []
klof_list, klod_list= [], []
for study in study_list:
    df = pd.read_pickle(path2octupoles_matching_data+study)
    for i in range(len(df.axy)):
        if i == 0:
            if study == study_list[0]: # for the case of negative ayy remove the first (ayy=0) entry in order not to be duplicated
                continue
            else:
                axy_list.append(df.axy[i])
                ayy_list.append(df.ayy[i])
                klof_list.append(df.klof[i])
                klod_list.append(df.klod[i])
        else:
            #if study == study_list[0]:
            #    df.axy, df.ayy, df.klof, df.klod = df.axy[::-1], df.ayy[::-1], df.klof[::-1], df.klod[::-1]
            axy_list.append(df.axy[i][0])
            ayy_list.append(df.ayy[i][0])
            klof_list.append(df.klof[i][0])
            klod_list.append(df.klod[i][0])
            if study == study_list[0]:
                axy_list, ayy_list, klof_list, klod_list = axy_list[::-1], ayy_list[::-1], klof_list[::-1], klod_list[::-1]

O3_lof = cmpt_octupole_coefficient(np.array(klof_list), Brho)
O3_lod = cmpt_octupole_coefficient(np.array(klod_list), Brho)

I_lof = list(cmpt_current_klof(O3_lof))
I_lod = list(cmpt_current_klod(O3_lod))


test_list = [[ayy_list][0], [axy_list][0], [klof_list][0], [klod_list][0], [O3_lof][0], [O3_lod][0], [I_lof][0], [I_lod][0]]
index_list = ['ayy [1/m]', 'axy [1/m]', 'klof [1/m^4]', 'klod [1/m^4]', 'O3_lof [T/m^3]', 'O3_lod [T/m^3]', 'I_lof [A]', 'I_lod [A]']

df1 = pd.DataFrame(test_list, index=index_list).T

# save data frame to pickle
save2pickle=True
if save2pickle:
	df1.to_pickle("summary_QpxQpy1_b3b5b7_270GeV_ayyScan_axyNoConstraint_LOD.pkl")
# test prin
print(df1[df1.keys()[2]])
