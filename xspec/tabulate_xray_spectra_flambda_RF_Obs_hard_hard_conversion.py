"""
plots the AGN spectra from the torus model for different obscuration levels.

output:
/home/comparat/wwwDir/eRoMok/taurus/*.png
"""
import xspec
import numpy as n
import sys
import astropy.units as uu
import astropy.constants as cc 
from scipy.interpolate import interp1d
from scipy.integrate import quad

xspec.Xset.cosmo = "67.77 0. 0.692885"
PL=1.9
#xspec.Xset.cosmo = "50. 0. 0.5"
#xspec.Response()

#s1 = xspec.Spectrum("arf01_100nmAl_200nmPI_sdtq.fits")
#b1 = s1.background
#r1 = s1.response
#arfFileName = r1.arf

nh_vals = 10**n.arange(-2.1,4.11,0.1)#0.05)
z_vals = n.arange(0., 4.01, 0.1)

nh_val = 1000.# nh_vals[0]
redshift = 2. # z_vals[0]

ll_05 = (cc.c * cc.h / (0.5*uu.keV).to(uu.J)).to(uu.AA)
ll_2 = (cc.c * cc.h / (2.0*uu.keV).to(uu.J)).to(uu.AA)
ll_10 = (cc.c * cc.h / (10*uu.keV).to(uu.J)).to(uu.AA)


#fs1 = xspec.FakeitSettings(os.path.join(os.environ['DARKSIM_DIR'], 'model', "arf01_100nmAl_200nmPI_sdtq.fits"), exposure = 1500.0)
#fs1.background = "none"

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
import numpy as n
import os
#torus_model = os.path.join(os.environ['DARKSIM_DIR'], 'model', 'torus1006.fits')
#print(torus_model)

## reproduce the model from Aird 2015


def get_spec(
	nh_val = 1
	,PL=1.9 
	,redshift=0.
	,norm1 = 0.98
	,norm2 = 0.02
	,norm3 = 1.
	,rel_refl= 1.
	,incl = n.cos(30.*n.pi/180.)
	):
	m1 = xspec.Model("zwabs*(cabs*zpowerlw*zhighect+pexrav)+zpowerlw")
	m1.setPars(
	nh_val,     #2    2   zwabs      nH         10^22    1.00000      +/-  0.0          
	redshift,   #3    2   zwabs      Redshift            0.0          frozen
	nh_val,     #4    3   cabs       nH         10^22    1.00000      +/-  0.0          
	PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
	redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
	norm1,      #7    4   zpowerlw   norm                1.00000      +/-  0.0          
	50.,        #8    5   zhighect   cutoffE    keV      10.0000      +/-  0.0          
	200.,       #9    5   zhighect   foldE      keV      15.0000      +/-  0.0          
	redshift,   #10    5   zhighect   Redshift            0.0          frozen
	PL,         #16    8   pexrav     PhoIndex            2.00000      +/-  0.0          
	200.,       #17    8   pexrav     foldE      keV      100.000      +/-  0.0          
	rel_refl,   #18    8   pexrav     rel_refl            0.0          +/-  0.0          
	redshift,   #19    8   pexrav     Redshift            0.0          frozen
	1.,         #20    8   pexrav     abund               1.00000      frozen
	1.,         #21    8   pexrav     Fe_abund            1.00000      frozen
	incl,       #22    8   pexrav     cosIncl             0.450000     frozen
	norm3,       #23    8   pexrav     norm                1.00000      +/-  0.0          
	PL,         #11    6   zpowerlw   PhoIndex            1.00000      +/-  0.0          
	redshift,   #12    6   zpowerlw   Redshift            0.0          frozen
	norm2      #13    6   zpowerlw   norm                1.00000      +/-  0.0          
	)
	kevs = n.arange(0.1, 50, 0.1) #
	kevs = 10**n.arange(-1.,n.log10(50),0.01)
	fluxes = []
	nPh = []
	for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
		xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
		#print(xspec.Xset.cosmo)
		#xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
		fluxes.append( m1.flux[0]    )
		nPh.append(m1.flux[3] )
		x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]
	# x  keV
	# y/dE erg/(cm$^2$ s keV)
	# nP number of photons
	# energy interval
	# E = hc/lambda
	# x * y / dE = nu * f_nu = ll * f_ll
	ll = (cc.c * cc.h / (x*uu.keV).to(uu.J)).to(uu.AA)
	f_ll = x * y / (dE * ll ) # erg/(cm$^2$ s A)
	return ll, f_ll, x, y/dE, nP, dE


def get_ratio_z(ll_0, f_ll_0, redshift):
	itp_0 = interp1d(ll_0, f_ll_0)
	itp_z2 = interp1d(ll_0*(1+redshift), f_ll_0)
	F_0 = quad(itp_0, ll_10.value, ll_2.value)[0]
	F_z2 = quad(itp_z2, ll_10.value, ll_2.value)[0]
	return F_0, F_z2

def get_curve(redshift = 6.):
	ratios = n.zeros_like(nh_vals)
	zs  = n.zeros_like(nh_vals)
	for ii, nh_val in enumerate(nh_vals):
		ll_0, f_ll_0, x, fnu, nP, dE = get_spec( nh_val, redshift = 0. )
		F_0, F_z2 = get_ratio_z(ll_0, f_ll_0, redshift)
		ratios[ii] = F_z2/F_0
		zs[ii] = redshift
	return zs, nh_vals, ratios 

ZZ = n.zeros(( len(z_vals), len(nh_vals)))
RR = n.zeros(( len(z_vals), len(nh_vals)))
NH = n.zeros(( len(z_vals), len(nh_vals)))

for jj, zz in enumerate( z_vals ) :
	zi,nhi,rri = get_curve(redshift = zz)
	ZZ[jj] = zi
	RR[jj] = rri
	NH[jj] = nhi

n.savetxt(
	os.path.join(os.environ['GIT_VS'], 'data', 'xray_k_correction', 'hard-2-hard-z.txt'),
	n.transpose([ 
		n.hstack(( ZZ )),
		n.hstack(( NH )),
		n.hstack(( RR ))
	]) )

