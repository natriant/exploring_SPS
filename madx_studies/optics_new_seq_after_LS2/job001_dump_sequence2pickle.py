from cpymad.madx import Madx
import pickle

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

twtable = mad.twiss()
#print(twtable.keys())

with open('twiss_s.pkl', 'wb') as fid:
    pickle.dump(twtable.s, fid)

with open('twiss_betx.pkl', 'wb') as fid:
    pickle.dump(twtable.betx, fid)

with open('twiss_bety.pkl', 'wb') as fid:
    pickle.dump(twtable.bety, fid)

with open('twiss_Dx.pkl', 'wb') as fid:
    pickle.dump(twtable.dx, fid)

with open('twiss_l.pkl', 'wb') as fid:
    pickle.dump(twtable.l, fid)
