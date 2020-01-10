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

nh_vals = 10**n.arange(-2,4+0.01,0.5)#0.05)
z_vals = n.arange(0.,6.1,0.5)
#z_vals = n.arange(0.,5.,0.1)
#10**n.arange(-3,0.76+0.01,0.25)#,0.025)

nh_val = 1
PL=1.9 
redshift=0.
f_scatter=0.02
norm1 = 1-f_scatter
norm2 = f_scatter
norm3 = 1.
rel_refl= -1.
incl = n.cos(30.*n.pi/180.)
norm_PR = 1.

def get_fraction_obs_RF_RF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 2.
	kev_max_erosita_RF = 10.
	m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val,   #1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		nh_val,   #2    2   plcabs     nH         10^22    1.00000      +/-  0.0          
		3.,       #3    2   plcabs     nmax       (scale)  1.00000      
		1.,       #4    2   plcabs     FeAbun              1.00000      frozen
		7.11,     #5    2   plcabs     FeKedge    KeV      7.11000      frozen
		PL,       #6    2   plcabs     PhoIndex            2.00000      +/-  0.0          
		50.,      #7    2   plcabs     HighECut   keV      95.0000      frozen
		200.,     #8    2   plcabs     foldE               100.000      frozen
		1.0,      #9    2   plcabs     acrit               1.00000      frozen
		0.0,      #10    2   plcabs     FAST       (scale)  0.0          
		redshift, #11    2   plcabs     Redshift            0.0          frozen
		norm1,    #12    2   plcabs     norm                1.00000      +/-  0.0          
		PL,       #13    3   pexrav     PhoIndex            2.00000      +/-  0.0          
		200.,     #14    3   pexrav     foldE      keV      100.000      +/-  0.0          
		rel_refl, #15    3   pexrav     rel_refl            0.0          +/-  0.0          
		redshift, #16    3   pexrav     Redshift            0.0          frozen
		1.,       #17    3   pexrav     abund               1.00000      frozen
		1.,       #18    3   pexrav     Fe_abund            1.00000      frozen
		incl,     #19    3   pexrav     cosIncl             0.450000     frozen
		norm_PR,  #20    3   pexrav     norm                1.00000      +/-  0.0          
		PL,       #21    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
		redshift, #22    4   zpowerlw   Redshift            0.0          frozen
		norm2, #23    4   zpowerlw   norm                1.00000      +/-  0.0          
		0.01 #24    5   TBabs      nH         10^22    1.00000      +/-  0.0          
		)
	m1.pexrav.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.TBabs.nH = 0.01
	m1.plcabs.nH = 0.01
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


def get_fraction_obsF_obsF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 0.4
	kev_max_erosita_RF = 2.4
	m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val,   #1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		nh_val,   #2    2   plcabs     nH         10^22    1.00000      +/-  0.0          
		3.,       #3    2   plcabs     nmax       (scale)  1.00000      
		1.,       #4    2   plcabs     FeAbun              1.00000      frozen
		7.11,     #5    2   plcabs     FeKedge    KeV      7.11000      frozen
		PL,       #6    2   plcabs     PhoIndex            2.00000      +/-  0.0          
		50.,      #7    2   plcabs     HighECut   keV      95.0000      frozen
		200.,     #8    2   plcabs     foldE               100.000      frozen
		1.0,      #9    2   plcabs     acrit               1.00000      frozen
		0.0,      #10    2   plcabs     FAST       (scale)  0.0          
		redshift, #11    2   plcabs     Redshift            0.0          frozen
		norm1,    #12    2   plcabs     norm                1.00000      +/-  0.0          
		PL,       #13    3   pexrav     PhoIndex            2.00000      +/-  0.0          
		200.,     #14    3   pexrav     foldE      keV      100.000      +/-  0.0          
		rel_refl, #15    3   pexrav     rel_refl            0.0          +/-  0.0          
		redshift, #16    3   pexrav     Redshift            0.0          frozen
		1.,       #17    3   pexrav     abund               1.00000      frozen
		1.,       #18    3   pexrav     Fe_abund            1.00000      frozen
		incl,     #19    3   pexrav     cosIncl             0.450000     frozen
		norm_PR,  #20    3   pexrav     norm                1.00000      +/-  0.0          
		PL,       #21    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
		redshift, #22    4   zpowerlw   Redshift            0.0          frozen
		norm2, #23    4   zpowerlw   norm                1.00000      +/-  0.0          
		0.01 #24    5   TBabs      nH         10^22    1.00000      +/-  0.0          
		)
	m1.pexrav.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.TBabs.nH = 0.01
	m1.plcabs.nH = 0.01
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
	m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val,   #1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		nh_val,   #2    2   plcabs     nH         10^22    1.00000      +/-  0.0          
		3.,       #3    2   plcabs     nmax       (scale)  1.00000      
		1.,       #4    2   plcabs     FeAbun              1.00000      frozen
		7.11,     #5    2   plcabs     FeKedge    KeV      7.11000      frozen
		PL,       #6    2   plcabs     PhoIndex            2.00000      +/-  0.0          
		50.,      #7    2   plcabs     HighECut   keV      95.0000      frozen
		200.,     #8    2   plcabs     foldE               100.000      frozen
		1.0,      #9    2   plcabs     acrit               1.00000      frozen
		0.0,      #10    2   plcabs     FAST       (scale)  0.0          
		redshift, #11    2   plcabs     Redshift            0.0          frozen
		norm1,    #12    2   plcabs     norm                1.00000      +/-  0.0          
		PL,       #13    3   pexrav     PhoIndex            2.00000      +/-  0.0          
		200.,     #14    3   pexrav     foldE      keV      100.000      +/-  0.0          
		rel_refl, #15    3   pexrav     rel_refl            0.0          +/-  0.0          
		redshift, #16    3   pexrav     Redshift            0.0          frozen
		1.,       #17    3   pexrav     abund               1.00000      frozen
		1.,       #18    3   pexrav     Fe_abund            1.00000      frozen
		incl,     #19    3   pexrav     cosIncl             0.450000     frozen
		norm_PR,  #20    3   pexrav     norm                1.00000      +/-  0.0          
		PL,       #21    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
		redshift, #22    4   zpowerlw   Redshift            0.0          frozen
		norm2, #23    4   zpowerlw   norm                1.00000      +/-  0.0          
		0.01 #24    5   TBabs      nH         10^22    1.00000      +/-  0.0          
		)
	m1.pexrav.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.TBabs.nH = 0.01
	m1.plcabs.nH = 0.01
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


def get_fraction_hard_obs_RF_ObsF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 2.0
	kev_max_erosita = 10.0
	kev_min_erosita_RF = 2.0/(1+redshift)
	kev_max_erosita_RF = 10.0/(1+redshift)
	m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val,   #1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		nh_val,   #2    2   plcabs     nH         10^22    1.00000      +/-  0.0          
		3.,       #3    2   plcabs     nmax       (scale)  1.00000      
		1.,       #4    2   plcabs     FeAbun              1.00000      frozen
		7.11,     #5    2   plcabs     FeKedge    KeV      7.11000      frozen
		PL,       #6    2   plcabs     PhoIndex            2.00000      +/-  0.0          
		50.,      #7    2   plcabs     HighECut   keV      95.0000      frozen
		200.,     #8    2   plcabs     foldE               100.000      frozen
		1.0,      #9    2   plcabs     acrit               1.00000      frozen
		0.0,      #10    2   plcabs     FAST       (scale)  0.0          
		redshift, #11    2   plcabs     Redshift            0.0          frozen
		norm1,    #12    2   plcabs     norm                1.00000      +/-  0.0          
		PL,       #13    3   pexrav     PhoIndex            2.00000      +/-  0.0          
		200.,     #14    3   pexrav     foldE      keV      100.000      +/-  0.0          
		rel_refl, #15    3   pexrav     rel_refl            0.0          +/-  0.0          
		redshift, #16    3   pexrav     Redshift            0.0          frozen
		1.,       #17    3   pexrav     abund               1.00000      frozen
		1.,       #18    3   pexrav     Fe_abund            1.00000      frozen
		incl,     #19    3   pexrav     cosIncl             0.450000     frozen
		norm_PR,  #20    3   pexrav     norm                1.00000      +/-  0.0          
		PL,       #21    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
		redshift, #22    4   zpowerlw   Redshift            0.0          frozen
		norm2, #23    4   zpowerlw   norm                1.00000      +/-  0.0          
		0.01 #24    5   TBabs      nH         10^22    1.00000      +/-  0.0          
		)
	m1.pexrav.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.TBabs.nH = 0.01
	m1.plcabs.nH = 0.01
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed


def get_fraction_hard_RF_soft_obsF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 2.0/(1+redshift)
	kev_max_erosita_RF = 10.0/(1+redshift)
	m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val,   #1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		nh_val,   #2    2   plcabs     nH         10^22    1.00000      +/-  0.0          
		3.,       #3    2   plcabs     nmax       (scale)  1.00000      
		1.,       #4    2   plcabs     FeAbun              1.00000      frozen
		7.11,     #5    2   plcabs     FeKedge    KeV      7.11000      frozen
		PL,       #6    2   plcabs     PhoIndex            2.00000      +/-  0.0          
		50.,      #7    2   plcabs     HighECut   keV      95.0000      frozen
		200.,     #8    2   plcabs     foldE               100.000      frozen
		1.0,      #9    2   plcabs     acrit               1.00000      frozen
		0.0,      #10    2   plcabs     FAST       (scale)  0.0          
		redshift, #11    2   plcabs     Redshift            0.0          frozen
		norm1,    #12    2   plcabs     norm                1.00000      +/-  0.0          
		PL,       #13    3   pexrav     PhoIndex            2.00000      +/-  0.0          
		200.,     #14    3   pexrav     foldE      keV      100.000      +/-  0.0          
		rel_refl, #15    3   pexrav     rel_refl            0.0          +/-  0.0          
		redshift, #16    3   pexrav     Redshift            0.0          frozen
		1.,       #17    3   pexrav     abund               1.00000      frozen
		1.,       #18    3   pexrav     Fe_abund            1.00000      frozen
		incl,     #19    3   pexrav     cosIncl             0.450000     frozen
		norm_PR,  #20    3   pexrav     norm                1.00000      +/-  0.0          
		PL,       #21    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
		redshift, #22    4   zpowerlw   Redshift            0.0          frozen
		norm2, #23    4   zpowerlw   norm                1.00000      +/-  0.0          
		0.01 #24    5   TBabs      nH         10^22    1.00000      +/-  0.0          
		)
	m1.pexrav.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
	xspec.AllModels.calcFlux(str(kev_min_erosita)+" "+str(kev_max_erosita))
	flux_obs = m1.flux[0]#/(kev_max_erosita-kev_min_erosita)
	# rest frame intrinsic flux
	xspec.AllModels.show()
	m1.TBabs.nH = 0.01
	m1.plcabs.nH = 0.01
	xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
	flux_intrinsic = m1.flux[0]#/(kev_max_erosita_RF-kev_min_erosita_RF)
	xspec.AllModels.show()
	fraction_observed = flux_obs / flux_intrinsic
	return fraction_observed

nh_val = 1
PL=1.9 
redshift=0.
f_scatter=0.02
norm1 = 1-f_scatter
norm2 = f_scatter
norm3 = 1.
rel_refl= -1.
incl = n.cos(30.*n.pi/180.)
f_scat_name = str(int(f_scatter*100)).zfill(3)

## rest-frame to observed frame conversion in the 0.5 to 2 keV band
#frac_RF_ObsF = n.array([n.array([get_fraction_obs_RF_ObsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
#z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
#nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
#n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_A15_RF_soft_ObsF_soft_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_ObsF))]), header='z log_nh fraction_observed')

## rest-frame to observed frame conversion in the 0.5 to 2 keV band
#frac_RF_RF = n.array([n.array([get_fraction_obs_RF_RF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
#z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
#nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
#n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_A15_RF_soft_RF_hard_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')

## rest-frame to observed frame conversion in the 2-10 keV band
frac_RF_RF = n.array([n.array([get_fraction_hard_obs_RF_ObsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_A15_RF_hard_Obs_hard_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')

# rest-frame 2-10 to observed frame 0.5-2
frac_RF_RF = n.array([n.array([get_fraction_hard_RF_soft_obsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
n.savetxt( os.path.join(os.environ['GIT_VS'],"data","xray_k_correction", "fraction_observed_A15_RF_hard_Obs_soft_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')
							



#nh_val,     #2    2   zwabs      nH         10^22    1.00000      +/-  0.0       
#redshift,   #3    2   zwabs      Redshift            0.0          frozen
#nh_val,     #4    3   cabs       nH         10^22    1.00000      +/-  0.0       
#PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0       
#redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
#norm1,      #7    4   zpowerlw   norm                1.00000      +/-  0.0       
#50.,        #8    5   zhighect   cutoffE    keV      10.0000      +/-  0.0       
#200.,       #9    5   zhighect   foldE      keV      15.0000      +/-  0.0       
#redshift,   #10    5   zhighect   Redshift            0.0          frozen
#PL,         #16    8   pexrav     PhoIndex            2.00000      +/-  0.0      
#200.,       #17    8   pexrav     foldE      keV      100.000      +/-  0.0      
#rel_refl,   #18    8   pexrav     rel_refl            0.0          +/-  0.0      
#redshift,   #19    8   pexrav     Redshift            0.0          frozen
#1.,         #20    8   pexrav     abund               1.00000      frozen
#1.,         #21    8   pexrav     Fe_abund            1.00000      frozen
#incl,       #22    8   pexrav     cosIncl             0.450000     frozen
#norm3,       #23    8   pexrav     norm                1.00000      +/-  0.0     
#PL,         #11    6   zpowerlw   PhoIndex            1.00000      +/-  0.0      
#redshift,   #12    6   zpowerlw   Redshift            0.0          frozen
#norm2      #13    6   zpowerlw   norm                1.00000      +/-  0.0       								