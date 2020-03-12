import os, sys, glob
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
from astropy.table import Table, Column
from astropy_healpix import healpy
import numpy as n

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as p

fig_dir = '/data/data/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn'

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmoMD = FlatLambdaCDM(
	H0=67.77 * u.km / u.s / u.Mpc,
	Om0=0.307115)  # , Ob0=0.048206)
h = 0.6777
L_box = 1000.0 / h
cosmo = cosmoMD

cluster_file = os.path.join(fig_dir, 'input_clusters.fits')
clu = Table.read(cluster_file)

prefix = 'solid-'
detection_file = os.path.join(fig_dir, 'master_IN_DET_solidlinks.fits')

prefix = 'primary-'
detection_file = os.path.join(fig_dir, 'master_IN_DET_primary.fits')

prefix = ''
detection_file = os.path.join(fig_dir, 'master_IN_DET.fits')
det_0 = Table.read(detection_file)

val, ids = n.unique(det_0['INPUT_XRAY_image_path'], return_index =True )
det = det_0[ids]

b10_C = ( abs( clu['g_lat'] ) > 10 ) 
b10_D = ( abs( det['INPUT_g_lat'] ) > 10 ) 




M500c_C = clu['CBP_EM0'][b10_C]
M500c_D = det['INPUT_CBP_EM0'][b10_D]

bins = n.arange(3, 8, 0.25)
p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
N_C = p.hist(M500c_C, bins=bins, histtype='step', label='simulated', lw=4)[0]
N_D = p.hist(M500c_D, bins=bins, histtype='step', label='detected', lw=2)[0]
p.xscale('log')
p.yscale('log')
p.xlim((bins[0], bins[-1]))
p.xlabel('-log10(EM(0))')
p.ylabel('Number of sources')
p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'EM0-hist.png'))
p.clf()

p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
p.plot(n.arange(3, 8, 0.25)[:-1] + 0.125, 1.*N_D/N_C )
frac = 1.*N_D/N_C
p.errorbar(x=n.arange(3, 8, 0.25)[:-1] + 0.125, y=frac, xerr=0.125, yerr=frac*N_D**(-0.5) )
#p.xscale('log')
#p.yscale('log')
#p.xlim((bins[0], bins[-1]))
p.ylim((0.3,1.01))
p.xlabel('-log10(EM(0))')
p.ylabel('Detected fraction')
#p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'EM0-hist-detected-fraction-png'))
p.clf()


M500c_C = clu['redshift_R'][b10_C]
M500c_D = det['INPUT_redshift_R'][b10_D]

bins = n.arange(0, 1.5, 0.05)
p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
N_C = p.hist(M500c_C, bins=bins, histtype='step', label='simulated', lw=4)[0]
N_D = p.hist(M500c_D, bins=bins, histtype='step', label='detected', lw=2)[0]
p.xscale('log')
p.yscale('log')
p.xlim((bins[0], bins[-1]))
p.xlabel('redshift')
p.ylabel('Number of sources')
p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'redshift_R-hist.png'))
p.clf()

p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
p.plot(n.arange(0, 1.5, 0.05)[:-1] + 0.025, 1.*N_D/N_C )
frac = 1.*N_D/N_C
p.errorbar(x=n.arange(0, 1.5, 0.05)[:-1] + 0.025, y=frac, xerr=0.025, yerr=frac*N_D**(-0.5) )
#p.xscale('log')
#p.yscale('log')
#p.xlim((bins[0], bins[-1]))
p.ylim((0.3, 1.01))
p.xlabel('redshift')
p.ylabel('Detected fraction')
#p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'redshift_R-hist-detected-fraction-png'))
p.clf()



FX_C = clu['CLUSTER_FX_soft'][b10_C]
FX_D = det['INPUT_CLUSTER_FX_soft'][b10_D]

bins = 10**n.arange(-15.25, -9, 0.25)
p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
N_C = p.hist(FX_C, bins=bins, histtype='step', label='simulated', lw=4)[0]
N_D = p.hist(FX_D, bins=bins, histtype='step', label='detected', lw=2)[0]
p.xscale('log')
p.yscale('log')
p.xlim((bins[0], bins[-1]))
p.xlabel('F_X')
p.ylabel('Number of sources')
p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'FX-hist.png'))
p.clf()

p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
p.plot(n.arange(-15.25, -9.25, 0.25) + 0.125, 1.*N_D/N_C )
frac = 1.*N_D/N_C
p.errorbar(x=n.arange(-15.25, -9.25, 0.25) + 0.125, y=frac, xerr=0.125, yerr=frac*N_D**(-0.5) )
#p.xscale('log')
#p.yscale('log')
#p.xlim((bins[0], bins[-1]))
p.ylim((0.3 ,1.01))
p.xlabel('F_X')
p.ylabel('Detected fraction')
#p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'FX-hist-detected-fraction-png'))
p.clf()



M500c_C = clu['HALO_M500c'][b10_C]
M500c_D = det['INPUT_HALO_M500c'][b10_D]

bins = 10**n.arange(13, 16, 0.25)
p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
N_C = p.hist(M500c_C, bins=bins, histtype='step', label='simulated', lw=4)[0]
N_D = p.hist(M500c_D, bins=bins, histtype='step', label='detected', lw=2)[0]
p.xscale('log')
p.yscale('log')
p.xlim((bins[0], bins[-1]))
p.xlabel('M_500c')
p.ylabel('Number of sources')
p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'M500c-hist.png'))
p.clf()

p.figure(0, (4,4))
p.axes([0.18, 0.16, 0.75, 0.75])
p.plot(n.arange(13, 16, 0.25)[:-1] + 0.125, 1.*N_D/N_C )
frac = 1.*N_D/N_C
p.errorbar(x=n.arange(13, 16, 0.25)[:-1] + 0.125, y=frac, xerr=0.125, yerr=frac*N_D**(-0.5) )
#p.xscale('log')
#p.yscale('log')
#p.xlim((bins[0], bins[-1]))
p.ylim((0.3 ,1.01))
p.xlabel('M_500c')
p.ylabel('Detected fraction')
#p.legend(loc=0)
p.grid()
p.savefig(os.path.join(fig_dir, prefix+'M500c-hist-detected-fraction-png'))
p.clf()
