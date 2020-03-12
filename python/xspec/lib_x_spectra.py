

"""
plots the AGN spectra from the torus model for different obscuration levels.

output:
/home/comparat/wwwDir/eRoMok/taurus/*.png
"""
#import xspec
import numpy as n
import sys
import astropy.units as uu
import astropy.constants as cc 
from scipy.interpolate import interp1d
from scipy.integrate import quad

E, sig_E = n.loadtxt(os.path.join(os.environ['GIT_VS'], 'data', 'xray_k_correction', 'sigma-e.txt'), unpack=True)

# E  : energy in keV
sigma_E = interp1d(E, sig_E)

# E  : energy in keV
# z  : redshift
# nH : H density in units of 1e22 cm-2
zwabs = lambda E, z, nH : n.e**(-nH * sigma_E(E) * (1+z) )

# E  : energy in keV
# z  : redshift
# a  : power-law slope
# K  : normalization
zpwerlaw = lambda E, z, PL, K : K*(E*(1+z))**(-PL)

# optically thin compton scattering, sigma t = 0.0079 almost not varying with E.
# E  : energy in keV
# nH : H density in units of 1e22 cm-2
cabs = lambda E, nH : n.e**(nH * 0.0079 )

spec = lambda E, z, nH, PL, K, K2 : zpwerlaw(E, z, PL, K)  * cabs(E, nH) *  zwabs(E, z, nH ) + zpwerlaw(E, z, PL, K2)  

def spectrum( E, nH_val ):
	return spec(E, z = 0.
	,nH = nH_val # 1e24
	,PL = 2.
	,K = 1.
	,K2 = .01*2.2  )

x_kevs = 10**n.arange(-1,2,0.01)

p.figure(1, (6,6))
y1 = spectrum(x_kevs, 0.01)
p.plot(x_kevs, y1, lw=2 , label='20')

y1 = spectrum(x_kevs, 1)
p.plot(x_kevs, y1, lw=2 , label='22')

y1 = spectrum(x_kevs, 100)
p.plot(x_kevs, y1, lw=2 , label='24')

y1 = spectrum(x_kevs, 10000)
p.plot(x_kevs, y1, lw=2 , label='26')


p.xscale('log')
p.yscale('log')
p.xlabel(r'E keV')
p.ylabel('$Flux$')
p.ylim((1e-6, 1e0))
p.legend(frameon=False, loc=0, fontsize=9)
p.tight_layout()
p.grid()
p.savefig('xSpectrum.png')
p.clf()




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


def get_PL(redshift=0., PL=2., norm1=1):
	m1 = xspec.Model("zpowerlw")
	m1.setPars(
	PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
	redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
	norm1     #7    4   zpowerlw   norm                1.00000      +/-  0.0             
	)
	kevs = 10**n.arange(-1.,n.log10(50),0.01)
	x_kevs = 10**(0.5 *(n.log10(kevs[:-1]) + n.log10(kevs[1:])))
	fluxes = []
	nPh = []
	for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
		xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
		fluxes.append( m1.flux[0]    )
		nPh.append(m1.flux[3] )
		x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]
	
	return x_kevs, y


def get_sigma_E(nh_val = 1,redshift=0., PL=2., norm1=1./(3.68922734e-11)):
	m1 = xspec.Model("zwabs*zpowerlw")
	m1.setPars(
	nh_val,     #2    2   zwabs      nH         10^22    1.00000      +/-  0.0          
	redshift,   #3    2   zwabs      Redshift            0.0          frozen
	PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
	redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
	norm1     #7    4   zpowerlw   norm                1.00000      +/-  0.0             
	)
	kevs = 10**n.arange(-1.,n.log10(50),0.01)
	x_kevs = 10**(0.5 *(n.log10(kevs[:-1]) + n.log10(kevs[1:])))
	fluxes = []
	nPh = []
	for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
		xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
		fluxes.append( m1.flux[0]    )
		nPh.append(m1.flux[3] )
		x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]
	
	sig_e = -n.log(y)
	return x_kevs, sig_e


x_kevs, y_pl = get_PL(norm1=1./(3.68922734e-11) )
x_kevs, sig_e = get_sigma_E( )
sel = (sig_e!=n.inf)
y_val = sig_e[sel]*x_kevs [sel]**3*100

itp = interp1d( 
	n.hstack((0.01, 0.05, 0.1, x_kevs[sel], 200 )),
	n.hstack((30., 40., 50., y_val ,  n.mean(y_val[-10:]) )) 
	)

kevs = 10**n.arange(-2.,2.1,0.01)

p.figure(1, (6,6))
p.plot(x_kevs, sig_e*x_kevs**3*100, lw=2 )
p.plot(kevs, itp(kevs))
p.xscale('log')
p.yscale('log')
p.xlabel(r'E keV')
p.ylabel('$\sigma_E E^3_{keV} [10^{-24} cm^{-2}]$')
#p.legend(frameon=False, loc=0, fontsize=9)
p.tight_layout()
p.grid()
p.savefig(os.path.join(os.environ['GIT_VS'], 'data', 'xray_k_correction', 'sigma-e.png'))
#p.savefig(os.path.join(os.environ['GIT_VS'], 'figures', 'MD10/agn/NH', 'xspec-zwabs-sigma-e-z0.png')
p.clf()

n.savetxt(os.path.join(os.environ['GIT_VS'], 'data', 'xray_k_correction', 'sigma-e.txt'), n.transpose([itp.x, itp.y*1e-24/(itp.x**3*100)]) )
