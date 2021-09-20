from cpymad.madx import Madx

mad = Madx()
#mad.chdir('./madx')
#mad.call('sps/beams/lhc_beam_injection.beamx')
mad.call('sps_thin_crabcavity.madx')

#call, file = 'sps/beams/lhc_beam_injection.beamx';

