'''
To select differnt options for the matching modify the following files accoridng to the point below: /madx/sps/cmd/sps_matching.cmd

1. axy constraint --> constraint, expr=axy=axy_val; (no comment out)
2. axy no contraint --> !constraint, expr=axy=axy_val; (comment out)
3. match for ayy > 0.0 --> constraint, expr=ayy=ayy_val;
4. match for ayy < 0.0 --> constraint, expr=ayy=-ayy_val;
5. match using both LOD and LOF --> 
	vary, name=KLOD, STEP=1.E-8;
	vary, name=KLOF, STEP=1.E-8;

6. match using only LOD -->
	vary, name=KLOD, STEP=1.E-8;
	!vary, name=KLOF, STEP=1.E-8;
'''

from cpymad.madx import Madx
import pandas as pd
import numpy as np


# MAD-X parameters dictionary
madx_settings = {'QH':26.13, 'QV':26.18, 'QPH':1.0, 'QPV':1.0}
seq_name = 'sps'
harmonic_number = 4620


ayy_list = np.arange(0, 22000.0, 2000)
print(ayy_list)

#ayy_list = [10000.0, 12000.0]

for my_ayy in ayy_list:
    mad = Madx()
    mad.options.echo = False
    mad.options.info = False
    mad.warn = False
    mad.chdir('./madx')
    mad.call('sps_thin_crabcavity.madx')

    for parameter in madx_settings:
        setting = madx_settings[parameter]
        mad.input(f'{parameter} = {setting};')

    mad.use(seq_name)

    
    # Include b3b5b7 in MBA and MBB
    #mad.call('./sps/cmd/sps_setMultipoles_upto7.cmd')
    #mad.input('exec, set_Multipoles_270GeV;')
    #mad.call('./sps/cmd/sps_assignMultipoles_upto7.cmd')
    #mad.input('exec, AssignMultipoles;')


    # Tune and Chromaticity matching
    mad.call('./sps/cmd/sps_matching.cmd')
    mad.input('exec, SPS_matchtunes(QH, QV);')
    mad.input('exec, SPS_setchroma_Q26(QPH, QPV);')
    mad.input('acta.31637, harmon=%d;'%harmonic_number)
    mad.input('exec, match_chroma(QPH ,QPV);')

    # Match the octupoles
    print('Start octupoles matching')
    ayy_val, axy_val = my_ayy, 0.0
    mad.input(f'exec, match_octupoles({ayy_val}, {axy_val});') # use this line, if the input is the ayy, axy coefficients and then matching to klof, klod


    table_list = list(mad.table) # list of existing table names
    print(table_list[-1])

    my_table = mad.table.matching_results # get table as dict like object
   
    my_dict = {}
    my_dict['klod'] = my_table.klod
    my_dict['klof'] = my_table.klof
    my_dict['axx'] = my_table.axx
    my_dict['ayy'] = my_table.ayy
    my_dict['axy'] = my_table.axy

    if my_ayy == ayy_list[0]:
        df = pd.DataFrame.from_dict(my_dict)
    else:
        df = df.append(my_dict, ignore_index=True)
print(df)
df.to_pickle('matching_results_QpxQpy1_nob3b5b7_200GeV_positive_ayy_lod.pkl')


