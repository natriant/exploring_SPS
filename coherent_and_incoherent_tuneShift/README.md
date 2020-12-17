1) ./madx/ --> containes the required files to run madx
2) ./tune_pickle_data/---> when the tune of a large distribution is computed the trackind is done in htcondor and only the computed tunes, Qx and Qy, are returned, and the used from the post analysis scripts.
	2b) SPSheadtail_CC_noise_cmptTune.py: script that computes the tune for large number of particles. To be submitted in htcondor.
