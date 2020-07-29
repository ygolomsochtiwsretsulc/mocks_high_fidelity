"""
Tabulates the obscured flux fraction using the torus model and dirrefent values of NH

output:
NH, redshift, obscured fraction

on ds54
pyCONDA

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

nh_val_galactic = 0.03 # 1e20 cm^-2 median value for abs(g_lat)>20 in H14PI maps.
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

def get_f_hard_RF_hard_obsF(nh_val, redshift=0, nh_val_galactic = 0.01):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 2.0*(1+redshift)
	kev_max_erosita = 10.0*(1+redshift)
	kev_min_erosita_RF = 2.0
	kev_max_erosita_RF = 10.0
	redshift=0.
	#
	#kev_min_erosita = 2.0
	#kev_max_erosita = 10.0
	#kev_min_erosita_RF = 2.0 / (1+redshift)
	#kev_max_erosita_RF = 10.0 / (1+redshift)
	#redshift=0.
	#
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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


def get_f_hard_obsF_soft_obsF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 2.0
	kev_max_erosita_RF = 10.0
	redshift=0.
	#m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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


def get_f_hard_obsF_soft_obsF_z0(nh_val):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	redshift=0
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 2.0
	kev_max_erosita_RF = 10.0
	#m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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


def get_f_hard_RF_soft_obsF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5*(1+redshift)
	kev_max_erosita = 2.0*(1+redshift)
	kev_min_erosita_RF = 2.0
	kev_max_erosita_RF = 10.0
	redshift=0.
	#m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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


## obs frame 2-10 to observed frame 0.5-2. redshift=0
frac_RF_RF = n.array([get_f_hard_obsF_soft_obsF_z0(nh_val) for nh_val in nh_vals]) 
n.savetxt( os.path.join(os.environ['GIT_AGN_MOCK'],"data","xray_k_correction", "v3_redshift0_fraction_observed_A15_RF_hard_Obs_soft_fscat_"+f_scat_name+".txt"), n.transpose([ 22 + n.log10(nh_vals), frac_RF_RF]), header='log_nh fraction_observed')

## rest-frame to observed frame conversion in the 2-10 keV band
frac_RF_RF = n.array([n.array([get_f_hard_RF_hard_obsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
n.savetxt( os.path.join(os.environ['GIT_AGN_MOCK'],"data","xray_k_correction", "v3_fraction_observed_A15_RF_hard_Obs_hard_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')

## rest-frame 2-10 to observed frame 0.5-2
frac_RF_RF = n.array([n.array([get_f_hard_obsF_soft_obsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
n.savetxt( os.path.join(os.environ['GIT_AGN_MOCK'],"data","xray_k_correction", "v3_fraction_observed_A15_Obs_hard_Obs_soft_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')

## obs frame 2-10 to observed frame 0.5-2
frac_RF_RF = n.array([n.array([get_f_hard_RF_soft_obsF(nh_val, redshift) for nh_val in nh_vals]) for redshift in z_vals ])
z_all = n.array([n.array([ redshift for nh_val in nh_vals]) for redshift in z_vals ])
nh_all = n.array([n.array([ nh_val for nh_val in nh_vals]) for redshift in z_vals ])
n.savetxt( os.path.join(os.environ['GIT_AGN_MOCK'],"data","xray_k_correction", "v3_fraction_observed_A15_RF_hard_Obs_soft_fscat_"+f_scat_name+".txt"), n.transpose([n.hstack((z_all)), 22 + n.log10(n.hstack((nh_all))), n.hstack((frac_RF_RF))]), header='z log_nh fraction_observed')


sys.exit()


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


def get_fraction_obs_RF_RF(nh_val, redshift=0):#, kev_min_erosita = 0.5, kev_max_erosita = 2.0):
	print(nh_val, redshift)
	kev_min_erosita = 0.5
	kev_max_erosita = 2.0
	kev_min_erosita_RF = 2.
	kev_max_erosita_RF = 10.
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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
	m1 = xspec.Model("TBabs(plcabs + zgauss + constant*powerlaw + pexrav*constant)")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val_galactic,   #   1    1   TBabs      nH         10^22    1.00000      +/-  0.0   ||||            0.01    -0.0001          0          0     100000      1e+06         
		nh_val,            #   2    2   plcabs     nH         10^22    1.00000      +/-  0.0   ||||            0.01       0.01      1e-06      1e-06     100000     100000         
		3.,                #   3    2   plcabs     nmax       (scale)  1.00000                 ||||               3
		1.,                #   4    2   plcabs     FeAbun              1.00000      frozen     ||||               1      -0.01          0          0         10         10
		7.11,              #   5    2   plcabs     FeKedge    KeV      7.11000      frozen     ||||            7.11    -0.0711          7          7         10         10
		PL,                #   6    2   plcabs     PhoIndex            2.00000      +/-  0.0   ||||             1.9     -0.019          0          0          3          3         
		95.,               #   7    2   plcabs     HighECut   keV      95.0000      frozen     ||||              95      -0.95       0.01          1        100        200
		300.,              #   8    2   plcabs     foldE               100.000      frozen     ||||             300         -3          1          1      1e+06      1e+06
		1.0,               #   9    2   plcabs     acrit               1.00000      frozen     ||||               1      -0.01          0          0          1          1
		0.0,               #  10    2   plcabs     FAST       (scale)  0.0                     ||||               0
		redshift,          #  11    2   plcabs     Redshift            0.0          frozen     ||||               0     -0.005     -0.999     -0.999         10         10
		1.,                #  12    2   plcabs     norm                1.00000      +/-  0.0   ||||               1      -0.01          0          0      1e+20      1e+24          
		6.4,               #  13    3   zgauss     LineE      keV      6.50000      +/-  0.0   ||||             6.4     -0.064          0          0      1e+06      1e+06          
		0.05,              #  14    3   zgauss     Sigma      keV      0.100000     +/-  0.0   ||||            0.05    -0.0005          0          0         10         20          
		redshift,          #  15    3   zgauss     Redshift            0.0          frozen     |||| = p11
		0.01,              #  16    3   zgauss     norm                1.00000      +/-  0.0   ||||            0.01     0.0001          0          0      1e+20      1e+24      
		0.02,              #  17    4   constant   factor              1.00000      +/-  0.0   ||||            0.02    -0.0002          0          0        0.1        0.1      
		PL,                #  18    5   powerlaw   PhoIndex            1.00000      +/-  0.0   |||| = p6      
		1.,                #  19    5   powerlaw   norm                1.00000      +/-  0.0   |||| = p12/(1. + p11)/(1./(1. + p11))^( - p6)      
		PL,                #  20    6   pexrav     PhoIndex            2.00000      +/-  0.0   |||| = p6          
		300,               #  21    6   pexrav     foldE      keV      100.000      +/-  0.0   ||||             300         -3          1          1      1e+06      1e+06          
		-1,                #  22    6   pexrav     rel_refl            0.0          +/-  0.0   ||||              -1      -0.01         -3         -3     -1e-09     -1e-09      
		redshift,          #  23    6   pexrav     Redshift            0.0          frozen     |||| = p11
		1.,                 #  24    6   pexrav     abund               1.00000      frozen     ||||               1      -0.01          0          0      1e+06      1e+06
		1.,					#  25    6   pexrav     Fe_abund            1.00000      frozen    ||||               1      -0.01          0          0      1e+06      1e+06
		0.45,				#  26    6   pexrav     cosIncl             0.450000     frozen    ||||            0.45    -0.0045       0.05       0.05       0.95       0.95
		1.,					#  27    6   pexrav     norm                1.00000      +/-  0.0  |||| = p19       
		1.					#  28    7   constant   factor              1.00000      +/-  0.0  ||||               1      -0.01       0.01       0.01        100        100       
		)
	#m1.pexrav.rel_refl.values
	m1.powerlaw.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
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


