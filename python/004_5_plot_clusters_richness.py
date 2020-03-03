"""
What it does
------------

Plots the cluster model, richness vs mass


"""
from astropy.table import Table, Column
from scipy.stats import scoreatpercentile
import glob
import sys
from astropy_healpix import healpy
import os
import time
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
from scipy.optimize import curve_fit
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import numpy as n
print('Create file with galaxies around clusters')
print('=> Abundance matching for magnitudes')
print('=> Red sequence colors')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


env =  sys.argv[1]
delta_crit = sys.argv[2]
print(env, delta_crit)
test_dir = os.path.join(os.environ[env])
all_catalogs = n.array(glob.glob(os.path.join(test_dir, 'fits', 'all_*_galaxiesAroundClusters.fit')))
all_aexp = n.array([ os.path.basename(el).split('_')[1] for el in all_catalogs ])
all_aexp.sort()

def mass_2_richness(M200c, redshift): 
	"""
	M200c to richness conversion using the scaling relation 
	Table 2, second line of Capasso et al. 2019 (1812.0609
	"""
	return 39.8 * (M200c/3e14)**(0.98) * ((1+redshift)/(1+0.18))**(-1.08)



fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'clusters', 'galaxies')
if os.path.isdir(fig_dir) == False:
	os.system('mkdir -p ' + fig_dir)

for aexp_str in all_aexp[::-1]:
	baseName = 'all_'+aexp_str # sys.argv[2]  # "all_0.62840"
	z_snap = 1./float(baseName.split('_')[1])-1.
	aexp_str = str(int(float(baseName.split('_')[1])*1e5)).zfill(6)
	print(env, baseName)
	test_dir = os.path.join(os.environ[env])

	# cosmology set up
	if env == "MD10" or env == "MD04":
		cosmoMD = FlatLambdaCDM(
			H0=67.77 * u.km / u.s / u.Mpc,
			Om0=0.307115)  # , Ob0=0.048206)
		h = 0.6777
		L_box = 1000.0 / h
		cosmo = cosmoMD
	if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
		cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
		h = 0.6774
		L_box = 1000.0 / h
		cosmo = cosmoUNIT

	# import all pathes
	path_2_CLU_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')
	path_2_CLU_SAT_catalog = os.path.join(test_dir, 'fits', baseName + '_galaxiesAroundClusters.fit')
	#path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT.fit')
	#path_2_CLU_SAT_RS_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT_RS.fit')
	# opens mock catalogs

	# galaxies
	#hdu_clu = fits.open(path_2_CLU_SAT_catalog)  # , mode='update')
	hdu_clu_i = Table.read(path_2_CLU_SAT_catalog)  # , mode='update')
	id_U = n.unique(hdu_clu_i['HOST_XRAY_image_path'], return_index = True )[1]
	hdu_clu = hdu_clu_i[id_U]
	zr_gal = hdu_clu['redshift_R']
	z_cluster = hdu_clu['HOST_redshift_R']
	DM = hdu_clu['mag_abs_r'] - hdu_clu['mag_r'] 
	mock_mag_abs = hdu_clu['mag_abs_r'] 
	#mock_mag_abs_sdss_u = hdu_clu['sdss_u'] + DM
	#mock_mag_abs_sdss_g = hdu_clu['sdss_g'] + DM
	#mock_mag_abs_sdss_r = hdu_clu['sdss_r'] + DM
	#mock_mag_abs_sdss_i = hdu_clu['sdss_i'] + DM
	#mock_mag_abs_sdss_z = hdu_clu['sdss_z'] + DM

	omega = lambda zz: cosmo.Om0*(1+zz)**3. / cosmo.efunc(zz)**2
	DeltaVir_bn98 = lambda zz : (18.*n.pi**2. + 82.*(omega(zz)-1)- 39.*(omega(zz)-1)**2.)/omega(zz)
	HOST_HALO_Mvir = hdu_clu['HOST_HALO_Mvir'] / h
	HOST_HALO_Rvir = hdu_clu['HOST_HALO_Rvir']
	HOST_HALO_M500c = hdu_clu['HOST_HALO_M500c'] / h
	HOST_HALO_R500c = (DeltaVir_bn98(z_snap)/500. * HOST_HALO_M500c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir
	HOST_HALO_M200c = hdu_clu['HOST_HALO_M200c'] / h
	HOST_HALO_R200c = (DeltaVir_bn98(z_snap)/200. * HOST_HALO_M200c / HOST_HALO_Mvir)**(1./3.)*HOST_HALO_Rvir

	zmin = n.min(zr_gal) 
	zmax = n.max(zr_gal) 
	print(zmin, '<z<', zmax)
	mean_redshift=(zmin+zmax)/2.
	def mass_2_richness_w_param(logM200c, M_slope=1, Norm=40): 
		"""
		M200c to richness conversion using the scaling relation 
		Table 2, second line of Capasso et al. 2019 (1812.0609
		"""
		return M_slope * logM200c + Norm
	ok = (hdu_clu['richness']>20)&(HOST_HALO_M200c>1e14)
	coeff = n.polyfit(n.log10(HOST_HALO_M200c[ok]), n.log10(hdu_clu['richness'][ok]), deg=1)
	mrange = n.arange(13, 16, 0.1)
	print(coeff)
	y_model = n.polyval(coeff, mrange)
	#popt, pcov = curve_fit(f = mass_2_richness_w_param, xdata = n.log10(HOST_HALO_M200c[ok]), ydata = n.log10(hdu_clu['richness'][ok]), p0=(1,40))#, sigma = hdu_clu['richness'][ok]**(0.5))
	#print(popt, pcov)
	fig_out = os.path.join(fig_dir, 'richness_mass_'+aexp_str+'.png')
	coeff_out = os.path.join(fig_dir, 'richness_mass_'+aexp_str+'.txt')
	n.savetxt(coeff_out, coeff)
	mdex = 0.1
	mbins = n.arange(-30, -10, mdex)
	p.figure(1, (6., 5.5))
	p.axes([0.15, 0.15, 0.8, 0.8])
	p.plot(HOST_HALO_M200c, hdu_clu['richness'], 'k,', label='mock', rasterized=True)
	p.plot(10**mrange, 10**y_model, 'r--', label='fit')
	##p.plot(mrange, mass_2_richness(mrange, z_snap), 'r--', label='C19')
	#p.plot(mrange, mass_2_richness_w_param(mrange, popt[0], popt[1]), 'b--', label='fit')

	p.title(r'$\bar{z}$=' + str(n.round(z_snap, 3)))
	p.xlabel(r'$M_{200c}$')
	p.ylabel('richness')
	p.yscale('log')
	p.xscale('log')
	p.legend(loc=0, frameon=False)
	p.grid()
	p.savefig(fig_out)
	p.clf()


