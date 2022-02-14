from cpymad.madx import Madx
import pandas as pd
import numpy as np


# MAD-X parameters dictionary
madx_settings = {'QH':26.13, 'QV':26.18, 'QPH':0.0, 'QPV':0.0}
seq_name = 'sps'
harmonic_number = 4620


ayy_list = np.arange(0, 21000.0, 1000)
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

    
    # Tune and Chromaticity matching
    mad.call('./sps/cmd/sps_matching.cmd')
    mad.input('exec, SPS_matchtunes(QH, QV);')
    mad.input('exec, SPS_setchroma_Q26(QPH, QPV);')
    mad.input('acta.31637, harmon=%d;'%harmonic_number)
    mad.input('exec, match_chroma(QPH ,QPV);')

    # Include b3b5b7 in MBA and MBB
    #mad.call('./sps/cmd/sps_setMultipoles_upto7.cmd')
    #mad.input('exec, set_Multipoles_270GeV;')
    #mad.call('./sps/cmd/sps_assignMultipoles_upto7.cmd')
    #mad.input('exec, AssignMultipoles;')

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
df.to_pickle('matching_results.pkl')


