"""
What it does
------------

Adds coordinates in degrees

Equatorial coordinates: Right Ascension, Declination
Galactic coordinates: latitude and longitude
Ecliptic coordinates: latitude and longitude

Redshift in real space and observed space
luminosity distance [cm] (in the cosmology of the simulation)

Foreground map values:
Planck extinction E(B-V) values
H14PI extinction H1 [cm-2] values

References
----------

Planck E(B-V) map. https://ui.adsabs.harvard.edu/abs/2014A&A...571A..11P
H14PI map. http://adsabs.harvard.edu/abs/2016A\%26A...594A.116H

command to run
--------------

python3 001_coordinates.py environmentVAR fileBasename

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

It will find the input file :
$environmentVAR/hlists/fits/${fileBasename}.fits

And will write outputs in h5 format:
$environmentVAR/hlists/fits/${fileBasename}_coordinates.fits

Figures (Optional)
if the variable 'make_figure' is set to True, then figures will be created in the git repo here :
$GIT_AGN_MOCK/figures/$environmentVAR/coordinates/
It makes relevant histograms and scatter plots for all columns present in the new file.

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, dustmaps, matplotlib


"""
# import python packages

import sys, os
print(int(sys.version[0]))
if int(sys.version[0])==2:
    import healpy
if sys.version[0]=='3':
    from astropy_healpix import healpy
from dustmaps.planck import PlanckQuery
from dustmaps.config import config
from astropy.coordinates import SkyCoord
from astropy import wcs
from scipy.interpolate import interp1d
import h5py
import astropy.io.fits as fits
from astropy.table import Table, Column
import numpy as n
import time
print('Adds coordinates')
t0 = time.time()

# import all pathes

env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "sat_0.62840"
print(env, baseName)
make_figure = True
#make_figure = False

# Data file dependencies
path_2_NH_map = '/data17s/darksim/observations/h1_maps/H14PI/asu.fit'

path_2_planck_EBV_map = '/data17s/darksim/observations/dust_maps/'
# this directory must contain this file a planck directory that contains :
# HFI_CompMap_ThermalDustModel_2048_R1.20.fits

test_dir = os.path.join(os.environ[env], 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')

# X-ray K-correction and attenuation curves
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


# Planck dust map setup
config['data_dir'] = path_2_planck_EBV_map
planck = PlanckQuery()
print("Python packages are loaded", time.time() - t0, 's')

# NH foreground map
NH_DATA2 = fits.open(path_2_NH_map)[1].data
nh_val = NH_DATA2['NHI']
print("NH map", time.time() - t0, 's')


print("interpolates z d_comoving and D_L", time.time() - t0, 's')
z_array = n.arange(0, 7.5, 0.001)
dc_to_z = interp1d(cosmo.comoving_distance(z_array), z_array)

d_L = cosmo.luminosity_distance(z_array)
dl_cm = (d_L.to(u.cm)).value
dL_interpolation = interp1d(z_array, dl_cm)


f1 = fits.open(path_2_light_cone)
xn = f1[1].data['x']
yn = f1[1].data['y']
zn = f1[1].data['z']
#print('x', xn.min(), xn.max())
#print('y', yn.min(), yn.max())
#print('z', zn.min(), zn.max())
vxn = f1[1].data['vx']
vyn = f1[1].data['vy']
vzn = f1[1].data['vz']
f1.close()

rr = (xn**2 + yn**2 + zn**2)**0.5
print('min rr', n.min(rr), time.time() - t0, 's')
print('max rr', n.max(rr), time.time() - t0, 's')

# angular coordinates
theta = n.arccos(zn / rr) * 180 / n.pi
phi = n.arctan2(yn, xn) * 180 / n.pi
ra = phi + 180.
dec = theta - 90.
print('ra', ra[:10], time.time() - t0, 's')
print('dec', dec[:10], time.time() - t0, 's')

# galactic and ecliptic coordinates
coords = SkyCoord(ra, dec, unit='deg', frame='icrs')
bb_gal = coords.galactic.b.value
ll_gal = coords.galactic.l.value
print('bb_gal', bb_gal[:10], time.time() - t0, 's')
print('ll_gal', ll_gal[:10], time.time() - t0, 's')

bb_ecl = coords.barycentrictrueecliptic.lat
ll_ecl = coords.barycentrictrueecliptic.lon
print('bb_ecl', bb_ecl[:10], time.time() - t0, 's')
print('ll_ecl', ll_ecl[:10], time.time() - t0, 's')

# line of sight, redshift
redshift_R = dc_to_z(rr)
vPara = (vxn * xn + vyn * yn + vzn * zn) / rr
rr_s = rr + vPara / cosmo.H(redshift_R).value
rr_s[rr_s <= 0] = rr[rr_s <= 0]
redshift_S = dc_to_z(rr_s)
print('redshift', redshift_R[:10], time.time() - t0, 's')

dL_cm = dL_interpolation(redshift_R)
print('dL_cm', dL_cm[:10], time.time() - t0, 's')

# H1 values
HEALPIX = healpy.ang2pix(
    1024,
    n.pi /
    2. -
    bb_gal *
    n.pi /
    180.,
    2 *
    n.pi -
    ll_gal *
    n.pi /
    180.)
NH = nh_val[HEALPIX]
print('NH', NH[:10], time.time() - t0, 's')

# E(B-V) values
ebv = planck(coords)
print('E(B-V)', ebv[:10], time.time() - t0, 's')

# Writes results
print('path_2_coordinate_file', path_2_coordinate_file)

t = Table()
t.add_column(Column(name='RA', data=ra, unit='deg'))
t.add_column(Column(name='DEC', data=dec, unit='deg'))

t.add_column(Column(name='g_lat', data=bb_gal, unit='deg'))
t.add_column(Column(name='g_lon', data=ll_gal, unit='deg'))

t.add_column(Column(name='ecl_lat', data=bb_ecl, unit='deg'))
t.add_column(Column(name='ecl_lon', data=ll_ecl, unit='deg'))

t.add_column(Column(name='redshift_R', data=redshift_R, unit=''))
t.add_column(Column(name='redshift_S', data=redshift_S, unit=''))

t.add_column(Column(name='dL', data=dL_cm, unit='cm'))
t.add_column(Column(name='nH', data=NH, unit='cm**(-2)'))
t.add_column(Column(name='ebv', data=ebv, unit='mag'))

t.write(path_2_coordinate_file, overwrite=True)
print('done', time.time() - t0, 's')


### NOW THE FIGURES ###
if make_figure:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': 14})
    import matplotlib.pyplot as p

    fig_dir = os.path.join(
        os.environ['GIT_AGN_MOCK'],
        'figures',
        env,
        'coordinates')
    if os.path.isdir(fig_dir) == False:
        os.system('mkdir -p ' + fig_dir)

    print('writes figures here: ', fig_dir)

    N_obj = len(redshift_R)
    rds = n.random.random(N_obj)
    sel = (rds < 1e6 / N_obj)

    # DEC vs redshift
    fig_out = os.path.join(fig_dir, 'Dec_redshift_R_' + baseName + '.png')
    XX = redshift_R[sel]
    YY = dec[sel]

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.plot(XX, YY, 'k,')
    p.title(baseName)
    p.xlabel('redshift')
    p.ylabel('Dec.')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    # DZ distribution
    fig_out = os.path.join(
        fig_dir,
        'redshift_R-minus-redshift_S_hist_' +
        baseName +
        '.png')

    dz = redshift_R - redshift_S

    X = n.log10(abs(dz[dz > 0]))

    p.figure(1, (5.5, 5.5))
    p.tight_layout()
    p.hist(X[X > -6], histtype='step', rasterized=True, lw=4, bins=100)
    p.title(baseName)
    p.xlabel('log10(|redshift R-redshiftS|)')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    #cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(E(B-V))$')
    #ax = p.gca()
    # ax.invert_xaxis()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'redshift_R_hist_' + baseName + '.png')
    p.figure(1, (5.5, 5.5))
    p.tight_layout()
    p.hist(redshift_R, histtype='step', rasterized=True, lw=4, bins=100)
    p.title(baseName)
    p.xlabel('redshift R')
    p.ylabel('Counts')
    p.yscale('log')
    p.grid()
    #cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(E(B-V))$')
    #ax = p.gca()
    # ax.invert_xaxis()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'redshift_S_hist_' + baseName + '.png')
    p.figure(1, (5.5, 5.5))
    p.tight_layout()
    p.hist(redshift_S, histtype='step', rasterized=True, lw=4, bins=100)
    p.title(baseName)
    p.xlabel('redshift S')
    p.ylabel('Counts')
    p.yscale('log')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'dL_hist_' + baseName + '.png')
    X = n.log10(dL_cm)
    p.figure(1, (5.5, 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4, bins=100)
    p.title(baseName)
    p.xlabel('dL [cm]')
    p.ylabel('Counts')
    p.yscale('log')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    # ebv plot

    fig_out = os.path.join(fig_dir, 'ra_dec_ebv_' + baseName + '.png')

    XX = ra[sel]
    YY = dec[sel]
    cc = n.log10(ebv[sel])

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.scatter(
        XX,
        YY,
        c=cc,
        s=3,
        marker='s',
        edgecolor='face',
        cmap='Paired',
        vmin=-3,
        vmax=2,
        rasterized=True)
    p.title(baseName)
    p.xlabel('R.A.')
    p.ylabel('Dec.')
    p.yscale('log')
    p.grid()
    cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(E(B-V))$')
    ax = p.gca()
    p.clf()

    fig_out = os.path.join(fig_dir, 'g_lat_g_lon_ebv_' + baseName + '.png')

    XX = ll_gal[sel]
    YY = bb_gal[sel]
    cc = n.log10(ebv[sel])

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.scatter(
        XX,
        YY,
        c=cc,
        s=3,
        marker='s',
        edgecolor='face',
        cmap='Paired',
        vmin=-3,
        vmax=2,
        rasterized=True)
    p.title(baseName)
    p.xlabel('galactic lon.')
    p.ylabel('galactic lat.')
    p.grid()
    cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(E(B-V))$')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'bb_ecl_ll_ecl_ebv_' + baseName + '.png')

    XX = ll_ecl[sel]
    YY = bb_ecl[sel]
    cc = n.log10(ebv[sel])

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.scatter(
        XX,
        YY,
        c=cc,
        s=3,
        marker='s',
        edgecolor='face',
        cmap='Paired',
        vmin=-3,
        vmax=2,
        rasterized=True)
    p.title(baseName)
    p.xlabel('ecliptic lon.')
    p.ylabel('ecliptic lat.')
    p.grid()
    cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(E(B-V))$')
    p.savefig(fig_out)
    p.clf()

    # NH plot

    fig_out = os.path.join(fig_dir, 'ra_dec_nh_' + baseName + '.png')

    XX = ra[sel]
    YY = dec[sel]
    cc = n.log10(NH[sel])

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.scatter(
        XX,
        YY,
        c=cc,
        s=3,
        marker='s',
        edgecolor='face',
        cmap='Paired',
        vmin=19.5,
        vmax=22,
        rasterized=True)
    p.title(baseName)
    p.xlabel('R.A.')
    p.ylabel('Dec.')
    p.grid()
    cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(n_H/cm^{-2})$')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'bb_gal_g_lon_nh_' + baseName + '.png')

    XX = ll_gal[sel]
    YY = bb_gal[sel]
    cc = n.log10(NH[sel])

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.scatter(
        XX,
        YY,
        c=cc,
        s=3,
        marker='s',
        edgecolor='face',
        cmap='Paired',
        vmin=19.5,
        vmax=22,
        rasterized=True)
    p.title(baseName)
    p.xlabel('galactic lon.')
    p.ylabel('galactic lat.')
    p.grid()
    cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(n_H/cm^{-2})$')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'bb_ecl_ll_ecl_nh_' + baseName + '.png')

    XX = ll_ecl[sel]
    YY = bb_ecl[sel]
    cc = n.log10(NH[sel])

    p.figure(0, (8., 5.5))
    p.tight_layout()
    p.scatter(
        XX,
        YY,
        c=cc,
        s=3,
        marker='s',
        edgecolor='face',
        cmap='Paired',
        vmin=19.5,
        vmax=22,
        rasterized=True)
    p.title(baseName)
    p.xlabel('ecliptic lon.')
    p.ylabel('ecliptic lat.')
    p.grid()
    cb = p.colorbar(shrink=0.7, label=r'$\log_{10}(n_H/cm^{-2})$')
    p.savefig(fig_out)
    p.clf()
