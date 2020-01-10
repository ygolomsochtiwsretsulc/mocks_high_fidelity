"""
What it does
------------

Computes the AGN model from Comparat et al. 2019.

It outputs X-ray and optical properties of AGNs

References
----------

 * Merloni et al. 2014 http://adsabs.harvard.edu/abs/2014MNRAS.437.3550M

 * Georgakakis et al. 2018 https://ui.adsabs.harvard.edu/\#abs/2018MNRAS.tmp.3272G

 * Comparat et al. 2019 https://ui.adsabs.harvard.edu/abs/2019MNRAS.tmp.1335C


Command to run
--------------

python3 003_0_agn.py environmentVAR fileBasename

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, extinction, matplotlib


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


env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "sat_0.62840"
print(env, baseName)
make_figure = True
make_figure = False


def get_a(baseName):
    alp = baseName.split('_')[1]
    print('a=', alp)
    return float(alp)


a_snap = get_a(baseName)

test_dir = os.path.join(os.environ[env], 'fits')

path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_agn_file = os.path.join(test_dir, baseName + '_agn.fits')

# link to X-ray K-correction and attenuation curves
path_2_hard_RF_obs_soft = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "xray_k_correction",
    "fraction_observed_A15_RF_hard_Obs_soft_fscat_002.txt")
path_2_RF_obs_hard = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    "data",
    "xray_k_correction",
    "fraction_observed_A15_RF_hard_Obs_hard_fscat_002.txt")
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
f2 = fits.open(path_2_coordinate_file)
zz = f2[1].data['redshift_R']
dL_cm = f2[1].data['dL']
galactic_NH = f2[1].data['nH']
galactic_ebv = f2[1].data['ebv']
N_galaxies = len(zz)
f2.close()

# computes the cosmological volume (full sky at that point)
zmin = n.min(zz)
zmax = n.max(zz)
z_mean = 0.5 * (zmin + zmax)
print(zmin, '<z<', zmax)
vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value)
DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
print('volume', vol, 'Mpc3')

print('opens galaxy file ', time.time() - t0)
f3 = fits.open(path_2_galaxy_file)
mass = f3[1].data['SMHMR_mass']  # log of the stellar mass
f3.close()

print('computes duty cycle ', time.time() - t0)
# duty cycle measured from Georgakakis 2017
f_duty = interp1d(
    n.array([0., 0.75, 2., 3.5, 10.1]),
    ##n.array([0.1, 0.2, 0.25, 0.25, 0.25])
    n.array([0.1, 0.2, 0.3, 0.3, 0.3])
)

if a_snap > 0.40:
    # all halos are written in the light cone at low redshift
    f_duty_realization = f_duty(zz)
if a_snap < 0.40:
    # only 30% of halos are written in the light cone at high redshifts
    f_duty_realization = f_duty(zz) / 0.3

active = (n.random.random(size=N_galaxies) <= f_duty_realization)
# ids to map to galaxy and halo files
ids_active = n.arange(N_galaxies)[active]

logm = mass[ids_active]
z = zz[ids_active]
dl_cm = dL_cm[ids_active]
n_agn = len(z)
#print(f_duty_realization, N_galaxies, n_agn, n_agn * 1. / N_galaxies)

# Hard LX Abundance Matching
# Equations 2 and 3 of Comparat et al. 2019
print('computes LX with HAM ', time.time() - t0)


def kz_h(z): return 10**(-4.03 - 0.19 * (1 + z))


def Ls_h(z): return 10**(44.84 - n.log10(((1 + 2.0) / (1 + z))
                                         ** 3.87 + ((1 + 2.0) / (1 + z))**(-2.12)))


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
fx_05_20 = 10**(lx_obs_frame_05_2) / (4 * n.pi * (dl_cm)**2.) / h**3
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
fx_2_10 = 10**(lx_obs_frame_2_10) / (4 * n.pi * (dl_cm)**2.) / h**3
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

attenuation = nh_law(galactic_NH[ids_active])
agn_rxay_flux_05_20_observed = fx_05_20 * attenuation
#print('agn_rxay_flux_05_20_observed',agn_rxay_flux_05_20_observed,time.time() - t0)


# optical extinction, Fitzpatrick 99
ebv_values = n.hstack((n.arange(0., 5., 0.01), 10**n.arange(1, 4, 0.1)))
ext_values = n.array([extinction.fitzpatrick99(
    n.array([6231.]), 3.1 * EBV, r_v=3.1, unit='aa')[0] for EBV in ebv_values])
ext_interp = interp1d(ebv_values, ext_values)
agn_rmag_observed = empirical_mag_r + ext_interp(galactic_ebv[ids_active])
#print('agn_rmag_observed', agn_rmag_observed, time.time() - t0)

# ===============================
# Writing results
# ===============================
print('writes', path_2_agn_file)
t = Table()
t.add_column(Column(name='ids_active', data=ids_active, unit=''))
t.add_column(Column(name='LX_hard', data=lx, unit='log10(L_X/[2-10keV, erg/s])'))
t.add_column(Column(name='LX_soft', data=lx_05_20, unit='log10(L_X/[0.5-2keV, erg/s])'))
t.add_column(Column(name='FX_soft', data=fx_05_20, unit='F_X / [0.5-2keV, erg/cm2/s]'))
t.add_column(Column(name='FX_soft_attenuated', data=agn_rxay_flux_05_20_observed, unit='F_X / [0.5-2keV, erg/cm2/s]'))
t.add_column(Column(name='FX_hard', data=fx_2_10, unit='F_X / [0.5-2keV, erg/cm2/s]'))
t.add_column(Column(name='logNH', data=logNH, unit='log10(nH/[cm-2])'))
t.add_column(Column(name='agn_type', data=opt_type, unit=''))
t.add_column(Column(name='random', data=random_number, unit=''))
t.add_column(Column(name='SDSS_r_AB', data=empirical_mag_r, unit='mag'))
t.add_column(Column(name='SDSS_r_AB_attenuated', data=agn_rmag_observed, unit='mag'))
t.write(path_2_agn_file, overwrite=True)
print('done', time.time() - t0, 's')


### OPTION: MAKE A SET OF FIGURES ###
if make_figure:
    print('makes figures ', time.time() - t0)

    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': 14})
    import matplotlib.pyplot as p

    fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'agn', )
    if os.path.isdir(fig_dir) == False:
        os.system('mkdir -p ' + fig_dir)

    h5_file = os.path.join(path_2_coordinate_file)
    f = h5py.File(h5_file, "r")
    zz = f['/coordinates/redshift_R']
    f.close()

    # AGN file
    h5_file = os.path.join(path_2_agn_file)
    f = h5py.File(h5_file, "r")
    AGN_ids_active = f['/AGN/ids_active']
    AGN_LX_hard = f['/AGN/LX_hard']
    AGN_LX_soft = f['/AGN/LX_soft']
    AGN_FX_soft = f['/AGN/FX_soft']
    AGN_FX_soft_attenuated = f['/AGN/FX_soft_attenuated']
    AGN_FX_hard = f['/AGN/FX_hard']
    AGN_logNH = f['/AGN/logNH']
    AGN_agn_type = f['/AGN/agn_type']
    AGN_random = f['/AGN/random']
    AGN_SDSS_r_AB = f['/AGN/SDSS_r_AB']
    AGN_SDSS_r_AB_attenuated = f['/AGN/SDSS_r_AB_attenuated']

    N_obj = len(AGN_ids_active)
    rds = n.random.random(N_obj)
    sel = (rds < 1e6 / N_obj)

    # 2D plots

    fig_out = os.path.join(fig_dir, 'AGN_LX_soft_vs_zz_' + baseName + '.png')

    X = zz[AGN_ids_active][sel]
    Y = AGN_LX_soft[sel]
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k+', rasterized=True)
    p.title(baseName)
    # p.yscale('log')
    p.xlabel('redshift R')
    p.ylabel('LX_soft')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'AGN_LX_hard_vs_mass_' + baseName + '.png')

    X = mass[AGN_ids_active][sel]
    Y = AGN_LX_hard[sel]
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k+', rasterized=True)
    p.title(baseName)
    # p.yscale('log')
    p.xlabel('mass')
    p.ylabel('LX_hard')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    # histograms

    fig_out = os.path.join(fig_dir, 'AGN_redshift_hist_' + baseName + '.png')
    X = zz[AGN_ids_active]

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=2, bins=100)
    p.title(baseName)
    p.xlabel('redshift')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(
        fig_dir,
        'AGN_SDSS_r_AB_attenuated_hist_' +
        baseName +
        '.png')
    X = AGN_SDSS_r_AB_attenuated

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(
        X,
        histtype='step',
        rasterized=True,
        lw=4,
        bins=n.arange(
            12,
            30,
            1.0))
    p.title(baseName)
    p.xlabel('AGN_SDSS_r_AB_attenuated')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'AGN_SDSS_r_AB_hist_' + baseName + '.png')
    X = AGN_SDSS_r_AB

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(
        X,
        histtype='step',
        rasterized=True,
        lw=4,
        bins=n.arange(
            12,
            30,
            1.0))
    p.title(baseName)
    p.xlabel('AGN_SDSS_r_AB')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'agn_type_hist_' + baseName + '.png')
    X = AGN_agn_type

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(
        X,
        histtype='step',
        rasterized=True,
        lw=4,
        bins=n.arange(
            10,
            23,
            1.))
    p.title(baseName)
    p.xlabel('agn_type')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'logNH_hist_' + baseName + '.png')
    X = AGN_logNH

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(
        X,
        histtype='step',
        rasterized=True,
        lw=4,
        bins=n.arange(
            20,
            26.5,
            0.5))
    p.title(baseName)
    p.xlabel('logNH')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'FX_hard_hist_' + baseName + '.png')
    X = n.log10(AGN_FX_hard)

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('FX_hard')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(
        fig_dir,
        'FX_soft_attenuated_hist_' +
        baseName +
        '.png')
    X = n.log10(AGN_FX_soft_attenuated)

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('FX_soft_attenuated')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'FX_soft_hist_' + baseName + '.png')
    X = n.log10(AGN_FX_soft)

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('FX_soft')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'LX_soft_hist_' + baseName + '.png')
    X = AGN_LX_soft

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('LX_soft')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'LX_hard_hist_' + baseName + '.png')

    X = AGN_LX_hard

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('LX_hard')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()
