from cpymad.madx import Madx

# MAD-X parameters dictionary
madx_settings = {'QH':26.13, 'QV':26.18, 'QPH':0.0, 'QPV':0.0}
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

# Include b3b5b7 in MBA and MBB
mad.call('./sps/cmd/sps_setMultipoles_upto7.cmd')
mad.input('exec, set_Multipoles_26GeV;')
mad.call('./sps/cmd/sps_assignMultipoles_upto7.cmd')
mad.input('exec, AssignMultipoles;')


# Tune matching. After the latest update the names of the quadrupoles changed. Therfore, to performe the tune matching correctly the following macro from the toolkit repository is needed.
# https://gitlab.cern.ch/acc-models/acc-models-sps/-/tree/2021/
mad.call('./sps/toolkit/macro.madx')
mad.input('exec, sps_match_tunes(QH,QV);')
#mad.call('./sps/cmd/sps_matching.cmd')
#mad.input('exec, SPS_matchtunes(QH, QV);')

# Chromaticity matching
mad.call('./sps/cmd/sps_matching.cmd')
mad.input('exec, SPS_setchroma_Q26(QPH, QPV);')
mad.input('acta.31637, harmon=%d;'%harmonic_number)
mad.input('exec, match_chroma(QPH ,QPV);')
print('here')



klod, klof = 0.0, 0.0
mad.input(f'klod = {klod};')
mad.input(f'klof = {klof};')
mad.call('./ptc/ptc.macro')
mad.input('exec, PTCchroma;')

print(f'klod {klod}, klof {klof}')

