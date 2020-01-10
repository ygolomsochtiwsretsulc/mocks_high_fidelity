"""
plots the AGN spectra from the torus model for different obscuration levels.

output:
/home/comparat/wwwDir/eRoMok/taurus/*.png
"""
import xspec
import numpy as n
import sys

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


p.figure(1, (7,4))

nh_val = 1
PL=1.9 
redshift=0.
f_scatter=0.02
norm1 = 1-f_scatter
norm2 = f_scatter
norm3 = 1.
rel_refl= 1.
incl = n.cos(30.*n.pi/180.)

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
kevs = 10**n.arange(-1.,n.log10(50),0.1)
fluxes = []
nPh = []
for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  #print(xspec.Xset.cosmo)
  #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  fluxes.append( m1.flux[0]    )
  nPh.append(m1.flux[3] )
x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

p.plot(x, y/dE, lw=4, label='all')


m1 = xspec.Model("zwabs*pexrav")

m1.setPars(
nh_val,     #2    2   zwabs      nH         10^22    1.00000      +/-  0.0          
redshift,   #3    2   zwabs      Redshift            0.0          frozen
PL,         #16    8   pexrav     PhoIndex            2.00000      +/-  0.0          
200.,       #17    8   pexrav     foldE      keV      100.000      +/-  0.0          
rel_refl,   #18    8   pexrav     rel_refl            0.0          +/-  0.0          
redshift,   #19    8   pexrav     Redshift            0.0          frozen
1.,         #20    8   pexrav     abund               1.00000      frozen
1.,         #21    8   pexrav     Fe_abund            1.00000      frozen
incl,       #22    8   pexrav     cosIncl             0.450000     frozen
norm3       #23    8   pexrav     norm                1.00000      +/-  0.0             
)


kevs = n.arange(0.1, 50, 0.1) #
kevs = 10**n.arange(-1.,n.log10(50),0.1)
fluxes = []
nPh = []
for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  #print(xspec.Xset.cosmo)
  #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  fluxes.append( m1.flux[0]    )
  nPh.append(m1.flux[3] )
x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

p.plot(x, y/dE, lw=2, ls='dashed', label='zwabs*pexrav')


m1 = xspec.Model("zwabs*cabs*zpowerlw*zhighect")

m1.setPars(
nh_val,     #2    2   zwabs      nH         10^22    1.00000      +/-  0.0          
redshift,   #3    2   zwabs      Redshift            0.0          frozen
nh_val,     #4    3   cabs       nH         10^22    1.00000      +/-  0.0          
PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
norm1,      #7    4   zpowerlw   norm                1.00000      +/-  0.0          
50.,        #8    5   zhighect   cutoffE    keV      10.0000      +/-  0.0          
200.,       #9    5   zhighect   foldE      keV      15.0000      +/-  0.0          
redshift   #10    5   zhighect   Redshift            0.0          frozen       
)


kevs = n.arange(0.1, 50, 0.1) #
kevs = 10**n.arange(-1.,n.log10(50),0.1)
fluxes = []
nPh = []
for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  #print(xspec.Xset.cosmo)
  #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  fluxes.append( m1.flux[0]    )
  nPh.append(m1.flux[3] )
x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

p.plot(x, y/dE, lw=1, ls='dashed', label='zwabs*cabs*zpowerlw*zhighect')

m1 = xspec.Model("zwabs*zpowerlw")

m1.setPars(
nh_val,     #2    2   zwabs      nH         10^22    1.00000      +/-  0.0          
redshift,   #3    2   zwabs      Redshift            0.0          frozen
PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
norm1     #7    4   zpowerlw   norm                1.00000      +/-  0.0             
)


kevs = n.arange(0.1, 50, 0.1) #
kevs = 10**n.arange(-1.,n.log10(50),0.1)
fluxes = []
nPh = []
for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  #print(xspec.Xset.cosmo)
  #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  fluxes.append( m1.flux[0]    )
  nPh.append(m1.flux[3] )
x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

p.plot(x, y/dE, lw=1, ls='dotted', label='zwabs*zpowerlw')

m1 = xspec.Model("zpowerlw")
m1.setPars(
PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
1.          
)

kevs = n.arange(0.1, 50, 0.1) #
kevs = 10**n.arange(-1.,n.log10(50),0.1)
fluxes = []
nPh = []
for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  #print(xspec.Xset.cosmo)
  #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  fluxes.append( m1.flux[0]    )
  nPh.append(m1.flux[3] )
x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

p.plot(x, y/dE, label='zpowerlw 1')

m1 = xspec.Model("zpowerlw")
m1.setPars(
PL,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
redshift,   #6    4   zpowerlw   Redshift            0.0          frozen
0.02          
)

kevs = n.arange(0.1, 50, 0.1) #
kevs = 10**n.arange(-1.,n.log10(50),0.1)
fluxes = []
nPh = []
for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
  xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  #print(xspec.Xset.cosmo)
  #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
  fluxes.append( m1.flux[0]    )
  nPh.append(m1.flux[3] )
x,y,nP,dE = (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

p.plot(x, y/dE, label='zpowerlw 0.02')

p.xlabel('keV observed frame')
p.xlim((0.2,50))
p.ylim((2e-12,1e-8))
p.yscale('log')
p.xscale('log')
p.title('redshift = '+str(redshift))
p.ylabel('erg/(cm$^2$ s keV)')
#p.ylabel('flux')
#p.ylabel('flux erg/cm2/s')
p.grid()
#p.axvline(0.5, label='eRosita observed band 0.5-2keV', color='m', ls='dashed')
#p.axvline(2, color='m', ls='dashed')
#p.axvline(2/(1+redshift), label='Rest-frame: 2-10keV band', color='k', ls='dotted')
#p.axvline(10./(1+redshift), color='k', ls='dotted')
p.legend(frameon=False, loc=0, fontsize=9)
p.tight_layout()
p.savefig(os.path.join(os.environ['GIT_VS'], 'figures','agn', 'torus', 'test.png'))
p.clf()

sys.exit()

m1 = xspec.Model("zpowerlw")
m1.setPars(
2.,         #5    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
0.,   #6    4   zpowerlw   Redshift            0.0          frozen
1.          
)

print(xspec.AllModels.calcFlux("0.5 2.5"))
print(xspec.AllModels.calcFlux("2 10"))


## model from Buchner et al. 2019  
torus_dir = os.path.join(os.environ['DARKSIM_DIR'], 'model', 'torus_buchner' )
torus_model = os.path.join(torus_dir, 'uxclumpy-cutoff.fits')
torus_model_omni = os.path.join(torus_dir, 'uxclumpy-cutoff-omni.fits')
torus_model_reflect = os.path.join(torus_dir, 'uxclumpy-cutoff-reflect.fits')
torus_model_transmit = os.path.join(torus_dir, 'uxclumpy-cutoff-transmit.fits')


def get_spectrum(nh_val, redshift):
  #m1 = xspec.Model("atable{"+torus_model+"} + zpowerlw")	
  #m1.setPars(nh_val, 2., 45., 87., redshift, 1., 2., redshift, 0.001)
  m1 = xspec.Model("atable{"+torus_model+"} + atable{"+torus_model_omni+"}")	
  m1.setPars(nh_val, PL, 400, 30, 0.3, 60., redshift, 1., nh_val, PL,  400, 30, 0.3, 60., redshift, 0.02)
  kevs = n.arange(0.1, 50, 0.1) #
  kevs = 10**n.arange(-1.,n.log10(50),0.1)
  fluxes = []
  nPh = []
  for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
    xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
    #print(xspec.Xset.cosmo)
    #xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
    fluxes.append( m1.flux[0]    )
    nPh.append(m1.flux[3] )
  return (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

def plot_spec(redshift):
  p.figure(1, (7,4))
  for nh_val in nh_vals:
    x,y,nP,dE=get_spectrum(nh_val, redshift)
    #p.plot(x*(1.+redshift), y, label=str(nh_val))
    #p.plot(x, nP/dE, label=str(nh_val))
    #print(x, y, nP, dE)
    p.plot(x, y/dE, label=str(n.round(n.log10(nh_val*10**22),1)))#"log(NH)="+

  p.xlabel('keV observed frame')
  p.xlim((0.2,20))
  p.yscale('log')
  p.xscale('log')
  p.title('redshift = '+str(redshift))
  p.ylabel('erg/(cm$^2$ s keV)')
  #p.ylabel('flux')
  #p.ylabel('flux erg/cm2/s')
  p.grid()
  #p.axvline(0.5, label='eRosita observed band 0.5-2keV', color='m', ls='dashed')
  #p.axvline(2, color='m', ls='dashed')
  #p.axvline(2/(1+redshift), label='Rest-frame: 2-10keV band', color='k', ls='dotted')
  #p.axvline(10./(1+redshift), color='k', ls='dotted')
  p.legend(frameon=False, loc=0, fontsize=9)
  p.tight_layout()
  p.savefig(os.path.join(os.environ['GIT_VS'], 'figures','agn', 'torus', 'spec_'+str(redshift)+'_omni_torus_scatter_2pc.png'))
  p.clf()

for zz in z_vals:
  plot_spec(zz)

sys.exit()

z, nh, fr = n.loadtxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_scatter01.txt"), unpack=True)

p.figure(2, (5,5))
p.scatter(z, nh, c=100.*fr, s=15, edgecolors='none', vmin=10, vmax=120)
cb = p.colorbar(shrink=0.7)
cb.set_label('fraction observed [%]')
p.xlabel('redshift')
#p.yscale('log')
#p.xscale('log')
#p.title('redshift = '+str(redshift))
p.ylabel('log(NH)')
#p.ylabel('flux')
#p.ylabel('flux erg/cm2/s')
p.grid()
#p.axvline(0.5, label='eRosita observed band 0.5-2keV', color='m', ls='dashed')
#p.axvline(2, color='m', ls='dashed')
#p.axvline(2/(1+redshift), label='Rest-frame: 2-10keV band', color='k', ls='dotted')
#p.axvline(10./(1+redshift), color='k', ls='dotted')
#p.legend(frameon=False, loc=0, fontsize=9)
p.savefig(os.path.join(os.environ['GIT_VS'], 'figures','agn', 'torus','fraction_observed_scatter01.png'))
p.clf()

#z, nh, fr_0001 = n.loadtxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_scatter0001.txt"), unpack=True)
#z, nh, fr_001 = n.loadtxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_scatter001.txt"), unpack=True)
#z, nh, fr_01 = n.loadtxt( os.path.join(os.environ['DARKSIM_DIR'], 'model', "fraction_observed_scatter01.txt"), unpack=True)


#p.figure(2, (5,5))
#p.scatter(z, nh, c=fr_001/fr_01, s=15, edgecolors='none', vmin=0, vmax=1)
#cb = p.colorbar(shrink=0.7)
#cb.set_label('f0.01 / f0.1')
#p.xlabel('redshift')
##p.yscale('log')
##p.xscale('log')
##p.title('redshift = '+str(redshift))
#p.ylabel('log(NH)')
##p.ylabel('flux')
##p.ylabel('flux erg/cm2/s')
#p.grid()
##p.axvline(0.5, label='eRosita observed band 0.5-2keV', color='m', ls='dashed')
##p.axvline(2, color='m', ls='dashed')
##p.axvline(2/(1+redshift), label='Rest-frame: 2-10keV band', color='k', ls='dotted')
##p.axvline(10./(1+redshift), color='k', ls='dotted')
##p.legend(frameon=False, loc=0, fontsize=9)
#p.savefig(os.path.join(os.environ['GIT_VS'], 'figures','agn', 'torus','fraction_observed_ratio_scatter01-001.png'))
#p.clf()

