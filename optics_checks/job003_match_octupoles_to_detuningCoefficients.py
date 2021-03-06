from cpymad.madx import Madx

# MAD-X parameters dictionary
madx_settings = {'QH':26.13, 'QV':26.18, 'QPH':1.0, 'QPV':1.0}
seq_name = 'sps'
harmonic_number = 4620

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


# Match the octupoles
print('Start octupoles matching')
ayy_val, axy_val = 1000, 0.0
mad.input(f'exec, match_octupoles({ayy_val}, {axy_val});') # use this line, if the input is the ayy, axy coefficients and then matching to klof, klod
