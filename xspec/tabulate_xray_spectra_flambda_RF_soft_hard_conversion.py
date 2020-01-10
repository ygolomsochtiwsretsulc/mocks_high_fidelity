

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

nh_vals = 10**n.arange(-2,4+0.01,2.)#0.05)
z_vals = n.arange(0., 4.1, 1.)

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


	#kevs = n.arange(0.1, 50, 0.1) #
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


def get_ratio(ll, f_ll_20):
	itp = interp1d(ll, f_ll_20/n.median(f_ll_20))
	F2_10 = quad(itp, ll_10.value, ll_2.value)[0]
	F05_2 = quad(itp, ll_2.value, ll_05.value)[0]
	return F05_2/F2_10

nh_vals = 10**n.arange(-2.1,4.11,0.1)#0.05)
ratios = n.zeros_like(nh_vals)
for ii, nh_val in enumerate(nh_vals):
	ll, f_ll, x, fnu, nP, dE = get_spec( nh_val )
	n.savetxt(os.path.join('xspectrum-flambda-'+str(n.round(n.log10(nh_val),1))+'.txt'), n.transpose([ll, f_ll, x, fnu, nP, dE]) )
	ratios[ii] = get_ratio(ll, f_ll)

n.savetxt(os.path.join(os.environ['GIT_VS'], 'data', 'xray_k_correction', 'hard-2-soft-z0.txt'), n.transpose([n.log10(nh_vals*1e22), ratios]) )

p.figure(1, (6,6))
p.plot(nh_vals*1e22, ratios)
p.xscale('log')
p.yscale('log')
p.xlabel(r'$\log_{10}(n_H)$')
p.ylabel('$F_{0.5-2}/F_{2-10}$')
#p.legend(frameon=False, loc=0, fontsize=9)
p.tight_layout()  
p.savefig(os.path.join(os.environ['GIT_VS'], 'figures', 'MD10/agn/NH', 'hard-2-soft-z0.png'))
p.clf()


ll, f_ll_21, x, fnu, nP, dE = get_spec( nh_val = 0.1, PL=1.9 )
ll, f_ll_22, x, fnu, nP, dE = get_spec( nh_val = 1, PL=1.9 )
ll, f_ll_23, x, fnu, nP, dE = get_spec( nh_val = 10, PL=1.9 )
ll, f_ll_24, x, fnu, nP, dE = get_spec( nh_val = 100, PL=1.9 )
ll, f_ll_25, x, fnu, nP, dE = get_spec( nh_val = 1000, PL=1.9 )
ll, f_ll_26, x, fnu, nP, dE = get_spec( nh_val = 10000, PL=1.9 )


p.figure(1, (7,4))
p.axvline(ll_05.value, label = '0.5 keV')
p.axvline(ll_2.value, label = '2 keV', ls='dashed')
p.axvline(ll_10.value, label = '10 keV', ls='dotted')
p.plot(ll, f_ll_20, label='20')
p.plot(ll, f_ll_21, label='21')
p.plot(ll, f_ll_22, label='22')
p.plot(ll, f_ll_23, label='23')
p.plot(ll, f_ll_24, label='24')
p.plot(ll, f_ll_25, label='25')
p.plot(ll, f_ll_26, label='26')
p.xscale('log')
p.yscale('log')
p.xlabel('wavelength [Angstrom]')
p.ylabel('flux density $f_\lambda$ [erg/(cm$^2$ s A)]')
p.legend(frameon=False, loc=0, fontsize=9)
p.tight_layout()  
p.savefig('ll_f_ll.png')
p.clf()
