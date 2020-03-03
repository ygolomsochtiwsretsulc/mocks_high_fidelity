"""
What it does
------------

Computes the AGN model from Comparat et al. 2019.


python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_2.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_4.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_6.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_8.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_10.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_12.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_14.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_16.fits
python 003_0_agn_EFEDS_WTHETA_model.py /home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_20.fits

python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_2.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_4.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_6.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_8.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_10.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_12.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_14.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_16.fits
python 003_0_agn_EFEDS_WTHETA_model.py /data37s/simulation_1/MD/MD_1.0Gpc/EFEDS_agn_fsat_20.fits



"""
import sys
import os
import time
import extinction
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from astropy.table import Table, Column
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
#import h5py
import numpy as n
print('Creates AGN mock catalogue ')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


env = 'MD10'
path_2_coordinate_file = sys.argv[1] # '/home/comparat/data/MultiDark/MD_1.0Gpc/EFEDS_agn_fsat_20.fits' # 
new_run = False
# link to X-ray K-correction and attenuation curves
#path_2_hard_RF_obs_soft = os.path.join(
    #os.environ['GIT_AGN_MOCK'],
    #"data",
    #"xray_k_correction",
    #"fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt")
#path_2_RF_obs_hard = os.path.join(
    #os.environ['GIT_AGN_MOCK'],
    #"data",
    #"xray_k_correction",
    #"fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt")
#path_2_NH_attenuation = os.path.join(
    #os.environ['GIT_AGN_MOCK'],
    #"data",
    #"xray_k_correction",
    #'gal_nh_ratio_relation_newg16.dat')

# link to X-ray K-correction and attenuation curves
path_2_hard_RF_obs_soft = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "xray_k_correction",
    "v2_fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt")
path_2_RF_obs_hard = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "xray_k_correction",
    "v2_fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt")
path_2_NH_attenuation = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "xray_k_correction",
    'gal_nh_ratio_relation_newg16.dat')

# simulation setup
if env[:2] == "MD" : # env == "MD04" or env == "MD40" or env == "MD10" or env == "MD25"
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env[:4] == "UNIT" : # == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR" or env == "MD40" or env == "MD10":
    scatter_0 = 1.4
if env == "MD04":
    scatter_0 = 1.0

print('opens coordinate file ', time.time() - t0)
f1 = Table.read( path_2_coordinate_file )
cen = (f1['pid']==-1)
sat = (cen==False)
zz = f1['redshift_R']
dL_cm = f1['dL']
galactic_NH = f1['nH']
galactic_ebv = f1['ebv']
mass = f1['SMHMR_mass']  # log of the stellar mass
#high_z = (zz>1.5)
if new_run :
	f1.add_column(Column(name='LX_hard', data=n.zeros_like(zz), unit='log10(L_X/[2-10keV, erg/s])'))
	f1.add_column(Column(name='LX_soft', data=n.zeros_like(zz), unit='log10(L_X/[0.5-2keV, erg/s])'))
	f1.add_column(Column(name='FX_soft', data=n.zeros_like(zz), unit='F_X / [0.5-2keV, erg/cm2/s]'))
	f1.add_column(Column(name='FX_soft_attenuated', data=n.zeros_like(zz), unit='F_X / [0.5-2keV, erg/cm2/s]'))
	f1.add_column(Column(name='FX_hard', data=n.zeros_like(zz), unit='F_X / [0.5-2keV, erg/cm2/s]'))
	f1.add_column(Column(name='logNH', data=n.zeros_like(zz), unit='log10(nH/[cm-2])'))
	f1.add_column(Column(name='agn_type', data=n.zeros_like(zz), unit=''))
	f1.add_column(Column(name='random', data=n.zeros_like(zz), unit=''))
	f1.add_column(Column(name='SDSS_r_AB', data=n.zeros_like(zz), unit='mag'))
	f1.add_column(Column(name='SDSS_r_AB_attenuated', data=n.zeros_like(zz), unit='mag'))

N_agn = len(zz)
N_agn_cen = len(zz[cen])
N_agn_sat = len(zz[sat])

print('obtained N, cen, sat, f_sat', N_agn, N_agn_cen, N_agn_sat, N_agn_sat*1./N_agn_cen)
print('N AGN', N_agn)

# computes the cosmological volume
area = 1600 # deg2
DZ = 0.1
z_bins = n.arange(0., 2., DZ)
##
for z_bins_i in z_bins:
	zmin = z_bins_i
	zmax = z_bins_i + DZ
	z_sel = (zz>=zmin) & (zz<zmax)
	z_mean = 0.5 * (zmin + zmax)
	print(zmin, '<z<', zmax)
	vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) * area * n.pi / 129600.
	DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
	print('volume', vol, 'Mpc3')
	logm = mass[z_sel]
	z = zz[z_sel]
	dl_cm = dL_cm[z_sel]
	n_agn = len(z)
	print('N AGN=', n_agn)
	# Hard LX Abundance Matching
	# Equations 2 and 3 of Comparat et al. 2019
	def kz_h(z): return 10**(-4.03 - 0.19 * (1 + z))

	def Ls_h(z): return 10**(44.84 - n.log10(((1 + 2.0) / (1 + z))** 3.87 + ((1 + 2.0) / (1 + z))**(-2.12)))

	def phi_h(L, z): return kz_h(z) / ((L / Ls_h(z))**0.48 + (L / Ls_h(z))**2.27)

	def scatter_z(x): return scatter_0 - 2 * x / 30.

	lsar_Zbar = n.zeros(n_agn)
	scatter = scatter_z(z_mean)

	# LF in the mock, starting parameters
	dlogf = 0.05
	Lbin_min = 36
	fbins = n.arange(Lbin_min, 48, dlogf)
	xf = fbins[:-1] + dlogf / 2.

	# theoretical number counts and LF
	N_obs_th = 1.0 * phi_h(10**xf, z_mean * n.ones_like(xf)) * vol * dlogf

	t1 = time.time()
	# select bins with a number of AGN greater than 1 and smaller than 2x the
	# total number of agn, we want to simulate
	bin_selection = (N_obs_th >= 0.5) & (N_obs_th < n_agn * 2.)
	# draw LX luminosities uniformly in each LX bin, the bins (dlogf = 0.05)
	# are small enough for a uniform sampling
	X_luminosities = n.hstack((
		n.array([n.random.uniform(low=aa, high=bb, size=cc)
				for aa, bb, cc in
				zip(fbins[:-1][bin_selection], fbins[1:][bin_selection], N_obs_th[bin_selection].astype('int') + 1)
				])
	))
	X_luminosities_sorted = X_luminosities[n.argsort(X_luminosities)]
	# print(X_luminosities_sorted)
	# scatter, then order the masses
	rds = norm.rvs(loc=0, scale=scatter, size=len(logm))
	M_scatt = logm + rds
	ids_M_scatt = n.argsort(M_scatt)
	# output numbers
	lx = n.zeros_like(logm)
	lx[ids_M_scatt] = X_luminosities_sorted[-n_agn:]
	lsar = n.zeros_like(lx)
	lsar[ids_M_scatt] = X_luminosities_sorted[-n_agn:] - logm[ids_M_scatt]

	# possibility: adjust redshift effect in the shell
	# lx = n.log10(DL_mean_z**2 * 10**lx / dl_cm**2)

	t2 = time.time()
	#print('HAM for LX needs N seconds/N agn= ', (t2 - t1) / n_agn)

	#print('lx', lx[:10], time.time() - t0)
	#print('lsar', lsar[:10], time.time() - t0)

	# ===============================
	# Obscured fractions
	# ===============================
	# model from equations 4-11, 12-15 of Comparat et al. 2019
	lx0 = 43.2


	def lxz(z): return lx0 + erf(z) * 1.2


	width = 0.6


	def thick_fraction_z(z): return 0.30  # + erf(z)*0.1


	def thin_fraction_max(LXhard): return 0.9 * (41 / LXhard)**0.5


	def thin_fraction_z(z): return thick_fraction_z(z) + 0.01 + erf(z / 4.) * 0.4


	def fraction_ricci(LXhard, z): return thin_fraction_z(z) + (thin_fraction_max(
		LXhard) - thin_fraction_z(z)) * (0.5 + 0.5 * erf((-LXhard + lxz(z)) / width))


	# initializes logNH
	logNH = n.zeros(n_agn)

	# obscuration, after the equations above
	randomNH = n.random.rand(n_agn)

	# unobscured 20-22
	#frac_thin = fraction_ricci(lsar, z)
	frac_thin = fraction_ricci(lx, z)
	thinest = (randomNH >= frac_thin)

	# thick obscuration, 24-26
	thick = (randomNH < thick_fraction_z(z))
	#thick = (randomNH < thick_fraction)

	# obscured 22-24
	obscured = (thinest == False) & (thick == False)

	# assigns logNH values randomly :
	logNH[thick] = n.random.uniform(24, 26, len(logNH[thick]))
	logNH[obscured] = n.random.uniform(22, 24, len(logNH[obscured]))
	logNH[thinest] = n.random.uniform(20, 22, len(logNH[thinest]))

	print('=====================  AGN fractions and numbers vs NH values =================')
	print(n_agn,
		len(thick.nonzero()[0]) * 1. / n_agn,
		len(obscured.nonzero()[0]) * 1. / n_agn,
		len(thinest.nonzero()[0]) * 1. / n_agn)

	# ===============================
	# Assigns flux
	# ===============================

	# hard X-ray 2-10 keV rest-frame ==>> 0.5-2 obs frame
	obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = n.loadtxt(
		path_2_hard_RF_obs_soft, unpack=True)
	obscuration_itp_H_S = interp2d(
		obscuration_z_grid,
		obscuration_nh_grid,
		obscuration_fraction_obs_erosita)

	NHS = n.arange(20, 26 + 0.05, 0.4)
	percent_observed_itp = interp1d(
		n.hstack((20 - 0.1, NHS, 26 + 0.1)),
		n.hstack((
			obscuration_itp_H_S(z_mean, 20.)[0],
			n.array([obscuration_itp_H_S(z_i, logNH_i)[0] for z_i, logNH_i in zip(z_mean * n.ones_like(NHS), NHS)]),
			obscuration_itp_H_S(z_mean, 26.)[0])))

	percent_observed_H_S = percent_observed_itp(logNH)

	lx_obs_frame_05_2 = n.log10(10**lx * percent_observed_H_S)
	fx_05_20 = 10**(lx_obs_frame_05_2) / (4 * n.pi * (dl_cm)**2.) # / h**3
	lx_05_20 = lx_obs_frame_05_2
	#print('fx_05_20', fx_05_20, time.time() - t0)
	#print('lx_05_20', lx_05_20, time.time() - t0)

	# hard X-ray 2-10 keV rest-frame ==>> 2-10 obs frame
	obscuration_z_grid, obscuration_nh_grid, obscuration_fraction_obs_erosita = n.loadtxt(
		path_2_RF_obs_hard, unpack=True)
	obscuration_itp_H_H = interp2d(
		obscuration_z_grid,
		obscuration_nh_grid,
		obscuration_fraction_obs_erosita)

	percent_observed_itp = interp1d(
		n.hstack((20 - 0.1, NHS, 26 + 0.1)),
		n.hstack((
			obscuration_itp_H_H(z_mean, 20.)[0],
			n.array([obscuration_itp_H_H(z_i, logNH_i)[0] for z_i, logNH_i in zip(z_mean * n.ones_like(NHS), NHS)]),
			obscuration_itp_H_H(z_mean, 26.)[0])))
	percent_observed_H_H = percent_observed_itp(logNH)

	lx_obs_frame_2_10 = n.log10(10**lx * percent_observed_H_H)
	fx_2_10 = 10**(lx_obs_frame_2_10) / (4 * n.pi * (dl_cm)**2.) # / h**3
	#print('fx_2_10', fx_2_10, time.time() - t0)
	#print('lx_obs_frame_2_10', lx_obs_frame_2_10, time.time() - t0)

	# Adds type 11, 12, 21, 22
	# Follows Merloni et al. 2014
	# equation 16 of Comparat et al. 2019


	def fraction_22p21_merloni(lx): return (
		0.5 + 0.5 * erf((-lx + 44.) / 0.9)) * 0.69 + 0.26


	def compute_agn_type(z, lx, logNH, fbins=fbins, n_agn=n_agn):
		"""
		Assigns a type to an AGN population

		parameters:
		- z: redshift
		- lx: hard X-ray luminosity (log10)
		- logNH: nH value (log10)

		return: array of AGN types
		"""
		# boundary between the 22 and the 21 populations
		limit = fraction_22p21_merloni((fbins[1:] + fbins[:-1]) * 0.5)
		# selection per obscuration intensity
		nh_21 = (logNH <= 22.)
		nh_23 = (logNH > 22.)  # &(logNH<=26.)
		# initiate columns to compute
		opt_type = n.zeros(n_agn).astype('int')
		rd = n.random.rand(n_agn)
		# compute histograms of LX for different obscurations
		nall = n.histogram(lx, fbins)[0]       # all
		nth = n.histogram(lx[nh_23], fbins)[0]  # thin
		nun = n.histogram(lx[nh_21], fbins)[0]  # unobscured
		fr_thk = nth * 1. / nall  # fraction of obscured
		fr_un = nun * 1. / nall  # fraction of unobscured
		# first get the type 12: NH absorption but optically unobscured
		# to be chosen in obscured population
		n_per_bin_12 = (fr_thk - limit) * nall
		sel_12 = (n.ones(len(z)) == 0)
		for bin_low, bin_high, num_needed, nn_un in zip(
				fbins[:-1], fbins[1:], n_per_bin_12.astype('int'), nth):
			if num_needed > 0 and nn_un > 0:
				frac_needed = num_needed * 1. / nn_un
				sel_12 = (sel_12) | (
					(lx > bin_low) & (
						lx < bin_high) & (nh_23) & (
						rd < frac_needed))
		t_12 = (nh_23) & (sel_12)
		# second the types 21
		# to be chosen in nun
		n_per_bin_21 = (-fr_thk + limit) * nall
		sel_21 = (n.ones(len(z)) == 0)
		for bin_low, bin_high, num_needed, nn_un in zip(
				fbins[:-1], fbins[1:], n_per_bin_21.astype('int'), nun):
			if num_needed > 0 and nn_un > 0:
				frac_needed = num_needed * 1. / nn_un
				sel_21 = (sel_21) | (
					(lx > bin_low) & (
						lx < bin_high) & (nh_21) & (
						rd < frac_needed))
		t_21 = (nh_21) & (sel_21)
		# finally the types 11 and 22
		t_11 = (nh_21) & (t_21 == False)
		t_22 = (nh_23) & (t_12 == False)
		opt_type[t_22] = 22
		opt_type[t_12] = 12
		opt_type[t_11] = 11
		opt_type[t_21] = 21
		return opt_type


	opt_type = compute_agn_type(z, lx, logNH)
	#print('opt_type', opt_type, time.time() - t0)

	# observed r-band magnitude from X-ray


	def r_mean(log_FX0520): return -2. * log_FX0520 - 7.


	def scatter_t1(n_agn_int): return norm.rvs(loc=0.0, scale=1.0, size=n_agn_int)


	random_number = n.random.rand(n_agn)
	empirical_mag_r = r_mean(n.log10(fx_05_20)) + scatter_t1(int(n_agn))
	#print('empirical_mag_r', empirical_mag_r, time.time() - t0)


	# ===============================
	# EXTINCTION
	# ===============================
	# x ray extinction from our Galaxy
	NH_DATA = n.loadtxt(path_2_NH_attenuation, unpack=True)
	nh_law = interp1d(
		n.hstack(
			(-10.**25, 10**n.hstack(
				(10., NH_DATA[0], 25)))), n.hstack(
					(1., 1., 1. / NH_DATA[1], 0.00001)))

	attenuation = nh_law(galactic_NH[z_sel])
	agn_rxay_flux_05_20_observed = fx_05_20 * attenuation
	#print('agn_rxay_flux_05_20_observed',agn_rxay_flux_05_20_observed,time.time() - t0)


	# optical extinction, Fitzpatrick 99
	ebv_values = n.hstack((n.arange(0., 5., 0.01), 10**n.arange(1, 4, 0.1)))
	ext_values = n.array([extinction.fitzpatrick99(
		n.array([6231.]), 3.1 * EBV, r_v=3.1, unit='aa')[0] for EBV in ebv_values])
	ext_interp = interp1d(ebv_values, ext_values)
	agn_rmag_observed = empirical_mag_r + ext_interp(galactic_ebv[z_sel])
	#print('agn_rmag_observed', agn_rmag_observed, time.time() - t0)

	# ===============================
	# Writing results
	# ===============================
	f1['LX_hard'][z_sel] = lx      
	f1['LX_soft'][z_sel] = lx_05_20
	f1['FX_soft'][z_sel] = fx_05_20
	f1['FX_soft_attenuated'][z_sel] = agn_rxay_flux_05_20_observed
	f1['FX_hard'][z_sel] = fx_2_10
	f1['logNH'][z_sel] = logNH
	f1['agn_type'][z_sel] = opt_type
	f1['random'][z_sel] = random_number
	f1['SDSS_r_AB'][z_sel] = empirical_mag_r
	f1['SDSS_r_AB_attenuated'][z_sel] = agn_rmag_observed

f1.write(path_2_coordinate_file, overwrite=True)
print('done', time.time() - t0, 's')

