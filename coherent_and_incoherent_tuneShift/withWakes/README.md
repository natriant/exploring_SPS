With wakefields ~5e5 particles are needed to obtain accurate results. Thus the simulation is to heavy to run locally and is submitted to condor. The tune computation is also done in htcondor, using NAFFlib.


1) create the bunch (grid in x-y) of 700^2 particles, using the script ../SPSheadtail_incoherentTuneShift_analyticalOnly.ipynb.
2) copy the bunch file along with the simulation script SPSheadtail_CC_noise_randomSeed.py in afs and run in htcondor.
	- output: Qx_list.pkl and Qy_list.pkl. The files contain the Qx and Qy as obtained from the tbt data for each particle.
	- copy the pickle files here to be used from the script: plot_footprint.py
3) Post process is perfromed here in plot_footprint.py
The footprint includes only the particles for which Jx/Jx_max + Jy/Jy_max < 1. 
The actions correspond to the initial actions.
