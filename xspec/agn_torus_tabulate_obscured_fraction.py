"""
Tabulates the obscured flux fraction using the torus model and dirrefent values of NH

output:
NH, redshift, obscured fraction

"""
import xspec
import numpy as n
import sys
import os
xspec.Xset.cosmo = "67.77 0. 0.692885"

nh_vals = 10**n.arange(-2,4+0.01,0.25)#0.05)
z_vals = n.hstack((0.001, n.arange(0.1,5.,0.1)))
#z_vals = n.arange(0.,5.,0.1)
#10**n.arange(-3,0.76+0.01,0.25)#,0.025)
PL=1.9
nh_val = 1000.# nh_vals[0]
redshift = 0. # z_vals[0]
f_scatter = 0.02
#torus_model = os.path.join(os.environ['DARKSIM_DIR'], 'model', 'torus1006.fits')
torus_dir = os.path.join(os.environ['DARKSIM_DIR'], 'model', 'torus_buchner' )
torus_model = os.path.join(torus_dir, 'uxclumpy-cutoff.fits')
torus_model_omni = os.path.join(torus_dir, 'uxclumpy-cutoff-omni.fits')
torus_model_reflect = os.path.join(torus_dir, 'uxclumpy-cutoff-reflect.fits')
torus_model_transmit = os.path.join(torus_dir, 'uxclumpy-cutoff-transmit.fits')


def get_fraction_hard_obs_RF_ObsF(nh_val, redshift=0):
	print(nh_val, redshift)
	kev_min_erosita = 2.0
	kev_max_erosita = 10.0
	kev_min_erosita_RF = 2.0/(1+redshift)
	kev_max_erosita_RF = 10.0/(1+redshift)
	m1 = xspec.Model("atable{"+torus_model+"} + atable{"+torus_model_omni+"}")	
	m1.setPars(nh_val,    # torus NH value in 0.01, 1000
		PL,       # torus photon Index, 2
		400,      # torus Ecut off 400 keV
		30,       # torus torus opening angle sigma. 5deg = thin, 90deg = isotropic
		0.3,      # torus CTK coverage: 0.4 ???
		40.,      # torus viewing angle
		redshift, # torus redshift
		1.,       # torus norm
		nh_val,   # scatt NH value in 0.01, 1000 (same as for torus)
		PL,       # scatt photon Index, 2
		400,      # scatt Ecut off 400 keV
		30,       # scatt torus opening angle sigma. 5deg = thin, 90deg = isotropic
		0.3,      # scatt CTK coverage: 0.4 ???
		40.,      # scatt viewing angle
		redshift, # scatt redshift
		f_scatter)    # scatt norm
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.torus.nH = 0.01
	m1.scat.nH = 0.01
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


def get_fraction_hard_RF_soft_obsF(nh_val, redshift=0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 2.0/(1+redshift)
	kev_max_erosita_RF = 10.0/(1+redshift)
	m1 = xspec.Model("atable{"+torus_model+"} + atable{"+torus_model_omni+"}")	
	m1.setPars(nh_val,    # torus NH value in 0.01, 1000
		PL,       # torus photon Index, 2
		400,      # torus Ecut off 400 keV
		30,       # torus torus opening angle sigma. 5deg = thin, 90deg = isotropic
		0.3,      # torus CTK coverage: 0.4 ???
		40.,      # torus viewing angle
		redshift, # torus redshift
		1.,       # torus norm
		nh_val,   # scatt NH value in 0.01, 1000 (same as for torus)
		PL,       # scatt photon Index, 2
		400,      # scatt Ecut off 400 keV
		30,       # scatt torus opening angle sigma. 5deg = thin, 90deg = isotropic
		0.3,      # scatt CTK coverage: 0.4 ???
		40.,      # scatt viewing angle
		redshift, # scatt redshift
		f_scatter)    # scatt norm
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.torus.nH = 0.01
	m1.scat.nH = 0.01
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed

for f_scatter in n.array([0.01, 0.02, 0.05, 0.1]):
  nh_val = 1
  PL=1.9 
  redshift=0.
  f_scat_name = str(int(f_scatter*100)).zfill(3)

  # rest-frame to observed frame conversion in the 2-10 keV band
  frac_RF_RF = n.array([n.array([get_fraction_hard_obs_RF_ObsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
  z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
  nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
  n.savetxt( os.path.join(os.environ['GIT_VS'], "data", "xray_k_correction", "fraction_observed_B19_RF_hard_Obs_hard_fscat_" + f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')

  # rest-frame 2-10 to observed frame 0.5-2
  frac_RF_RF = n.array([n.array([get_fraction_hard_RF_soft_obsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
  z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
  nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
  n.savetxt( os.path.join(os.environ['GIT_VS'], "data", "xray_k_correction", "fraction_observed_B19_RF_hard_Obs_soft_fscat_" + f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')


sys.exit()


print(torus_model)
def get_fraction_obs_RF_RF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
  print(nh_val, redshift)
  kev_min_erosita = 0.5
  kev_max_erosita = 2.0
  kev_min_erosita_RF = 2.
  kev_max_erosita_RF = 10.
  m1 = xspec.Model("atable{"+torus_model+"} + atable{"+torus_model_omni+"}")	
  m1.setPars(nh_val,    # torus NH value in 0.01, 1000
	      PL,       # torus photon Index, 2
	      400,      # torus Ecut off 400 keV
	      30,       # torus torus opening angle sigma. 5deg = thin, 90deg = isotropic
	      0.3,      # torus CTK coverage: 0.4 ???
	      40.,      # torus viewing angle
	      redshift, # torus redshift
	      1.,       # torus norm
	      nh_val,   # scatt NH value in 0.01, 1000 (same as for torus)
	      PL,       # scatt photon Index, 2
	      400,      # scatt Ecut off 400 keV
	      30,       # scatt torus opening angle sigma. 5deg = thin, 90deg = isotropic
	      0.3,      # scatt CTK coverage: 0.4 ???
	      40.,      # scatt viewing angle
	      redshift, # scatt redshift
	      0.02)    # scatt norm
  xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
  flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
  # rest frame intrinsic flux
  xspec.AllModels.show()
  m1.torus.nH = 0.01
  m1.scat.nH = 0.01
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
  xspec.AllModels.show()
  fraction_observed = flux_obs / flux_intrinsic
  return fraction_observed


def get_fraction_obs_RF_ObsF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
  print(nh_val, redshift)
  kev_min_erosita = 0.5
  kev_max_erosita = 2.0
  kev_min_erosita_RF = 0.5/(1+redshift)
  kev_max_erosita_RF = 2./(1+redshift)
  m1 = xspec.Model("atable{"+torus_model+"} + atable{"+torus_model_omni+"}")	
  m1.setPars(nh_val,    # torus NH value in 0.01, 1000
	      PL,       # torus photon Index, 2
	      400,      # torus Ecut off 400 keV
	      30,       # torus torus opening angle sigma. 5deg = thin, 90deg = isotropic
	      0.3,      # torus CTK coverage: 0.4 ???
	      40.,      # torus viewing angle
	      redshift, # torus redshift
	      1.,       # torus norm
	      nh_val,   # scatt NH value in 0.01, 1000 (same as for torus)
	      PL,       # scatt photon Index, 2
	      400,      # scatt Ecut off 400 keV
	      30,       # scatt torus opening angle sigma. 5deg = thin, 90deg = isotropic
	      0.3,      # scatt CTK coverage: 0.4 ???
	      40.,      # scatt viewing angle
	      redshift, # scatt redshift
	      0.02)    # scatt norm
  xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
  flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
  # rest frame intrinsic flux
  xspec.AllModels.show()
  m1.torus.nH = 0.01
  m1.scat.nH = 0.01
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
  xspec.AllModels.show()
  fraction_observed = flux_obs / flux_intrinsic
  return fraction_observed

## rest-frame to observed frame conversion in the 0.5 to 2 keV band
#frac_RF_ObsF = n.array([n.array([get_fraction_obs_RF_ObsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
#z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
#nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
##n.savetxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_scatter001.txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac))]), header='z log_nh fraction_observed')
#n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_scatter_1pc_RF_soft_ObsF_soft.txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_ObsF))]), header='z log_nh fraction_observed')

# rest-frame conversion from 2 to 10 keV band to the 0.5 to 2 keV band
frac_RF_RF = n.array([n.array([get_fraction_obs_RF_RF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
#n.savetxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_scatter001.txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac))]), header='z log_nh fraction_observed')
n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_scatter_1pc_RF_soft_RF_hard.txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')
#n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_scatter_1pc_RF_soft_RF_hard.txt"), n.transpose([22 + n.log10(nh_vals), frac_RF_RF]), header='log_nh fraction_observed')

sys.exit()


# PREVIOUS TORUS MODEL

#print(torus_model)
#def get_fraction_obs(nh_val, redshift, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
    #print(nh_val, redshift)
    #kev_min_erosita = 0.5#*(1+redshift)
    #kev_max_erosita = 2.0#*(1+redshift)
    #kev_min_erosita_RF = 2./(1+redshift)
    #kev_max_erosita_RF = 10./(1+redshift)
    #m1 = xspec.Model("atable{"+torus_model+"} + zpowerlw")	
    #m1.setPars(nh_val, PL, 45., 87., redshift, 1., PL, redshift, 0.01)
    ##m1.torus.nH = nh_val
    ##m1.torus.Theta_tor = 45.
    ##m1.torus.Theta_inc = 87.
    ##m1.torus.z = redshift
    ##m1.torus.norm = 1.
    ##m1.zpowerlw.PhoIndex=2.
    ##m1.zpowerlw.Redshift=redshift
    ##m1.zpowerlw.norm=0.01
    ### observed frame attenuated flux
    #xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
    #flux_obs = m1.flux[0]
    ## rest frame intrinsic flux
    #xspec.AllModels.show()
    #m1.torus.nH = 0.01
    #xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
    #flux_intrinsic = m1.flux[0]
    #xspec.AllModels.show()
    #fraction_observed = flux_obs / flux_intrinsic
    #return fraction_observed

#frac = n.array([n.array([get_fraction_obs(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])

#z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
#nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])

#n.savetxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_torus1006.txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac))]), header='z log_nh fraction_observed')

