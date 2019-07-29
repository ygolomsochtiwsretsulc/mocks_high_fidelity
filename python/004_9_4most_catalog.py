"""
Creates a fits catalog containing the 4FS input columns.

Create the ID array
remove unwanted column
check rulesets and templates names are correct

"""
import glob
import sys
from astropy_healpix import healpy
import os
from scipy.special import gammainc  # , gamma,  gammaincinv, gammaincc
from scipy.stats import scoreatpercentile
import pandas as pd  # external package
from scipy.special import erf
from astropy.coordinates import SkyCoord
import astropy.constants as cc
import astropy.io.fits as fits
from astropy.table import Table, Column
import astropy.units as u
import numpy as n
print('CREATES 4FS FITS FILES')
print('------------------------------------------------')
print('------------------------------------------------')

# astropy

#from lib_model_agn import *
#option = sys.argv[1]
#option = 'SNR3'
#option = 'eRASS8'
#option = 'eRASS3'
area_per_cat = healpy.nside2pixarea(8, degrees=True)

env = 'MD10'  # sys.argv[1]

# simulation setup
if env == "MD10" or env == "MD04":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

root_dir = os.path.join(os.environ[env])

dir_2_eRO_gal = os.path.join(root_dir, "cat_eRO_CLU_RS")
dir_2_GAL_all = os.path.join(root_dir, 'cat_GALAXY_all')
dir_2_GAL_sat = os.path.join(root_dir, 'cat_GALAXY_sat')

dir_2_4MOST = os.path.join(root_dir, "cat_CLU_4MOST")

if os.path.isdir(dir_2_4MOST) == False:
    os.system('mkdir -p ' + dir_2_4MOST)

#N_pixels = healpy.nside2npix(8)
# for HEALPIX_id in n.arange(N_pixels):
HEALPIX_id = 350
path_2_eRO_gal_catalog = os.path.join(
    dir_2_eRO_gal, str(HEALPIX_id).zfill(6) + '.fit')
path_2_GAL_all_catalog = os.path.join(
    dir_2_GAL_all, str(HEALPIX_id).zfill(6) + '.fit')
path_2_GAL_sat_catalog = os.path.join(
    dir_2_GAL_sat, str(HEALPIX_id).zfill(6) + '.fit')
path_2_4MOST_catalog = os.path.join(
    dir_2_4MOST,
    '4MOST_' +
    str(HEALPIX_id).zfill(6) +
    '.fit')

# to select BCG and red gal :
hd_clu = fits.open(path_2_eRO_gal_catalog)[1].data
# to select low z filament sample :
hd_all = fits.open(path_2_GAL_all_catalog)[1].data
hd_sat = fits.open(path_2_GAL_sat_catalog)[1].data

# =============================
# =============================
# eROSITA SAMPLE
# apply richness limit
# =============================
# =============================


area_ero = (abs(hd_clu['g_lat']) > 20) & (
    hd_clu['ecl_lon'] > 180) & (hd_clu['dec'] < 5)
r_10 = (hd_clu['richness'] > 40)  # & (area_ero)
bcg = (hd_clu['comoving_distance_to_cluster_in_rvir'] == 0.) & (
    hd_clu['HALO_M500c'] > 1e14)  # & (area_ero)
cgal = (r_10) & (
    hd_clu['galaxy_mag_r'] < 21.5) & (
        hd_clu['is_quiescent']) & (
            bcg == False)  # & (area_ero)

print(len(hd_clu['galaxy_mag_r'][bcg]) / area_per_cat)
print(len(hd_clu['galaxy_mag_r'][cgal]) / area_per_cat)

# =============================
# =============================
# LOWZ galaxy survey
# =============================
# =============================

area_all = (abs(hd_all['g_lat']) > 20) & (
    hd_all['ecl_lon'] > 180) & (hd_all['dec'] < 5)
z_sel_all = (hd_all['redshift_R'] < 0.3) & (
    hd_all['galaxy_stellar_mass'] > 10.8)  # & (area_all)

area_sat = (abs(hd_sat['g_lat']) > 20) & (
    hd_sat['ecl_lon'] > 180) & (hd_sat['dec'] < 5)
z_sel_sat = (hd_sat['redshift_R'] < 0.3) & (
    hd_sat['galaxy_stellar_mass'] > 10.8)  # & (area_sat)

print(len(hd_all['redshift_R'][z_sel_all]) / area_per_cat, 'cen per deg2')
print(len(hd_sat['redshift_R'][z_sel_sat]) / area_per_cat, 'sat per deg2')

# =============================
# =============================
# concatenate sub surveys
# =============================
# =============================

# 4 sub samples to concatenate


def get_qty(qty_name):
    return n.hstack(
        (hd_clu[qty_name][bcg],
         hd_clu[qty_name][cgal],
         hd_all[qty_name][z_sel_all],
         hd_sat[qty_name][z_sel_sat]))


subsurvey_id = n.hstack(
    (n.ones_like(
        hd_clu['ra'][bcg]),
        n.ones_like(
            hd_clu['ra'][cgal]) *
     2,
     n.ones_like(
        hd_all['ra'][z_sel_all]) *
     3,
     n.ones_like(
        hd_sat['ra'][z_sel_sat]) *
     3)).astype('int')

ra_array = get_qty('ra')
N_targets = len(ra_array)
dec_array = get_qty('dec')
galactic_ebv_array = get_qty('galactic_ebv')
z_array = get_qty('redshift_R')
stellar_mass = get_qty('galaxy_stellar_mass')

name_array = n.array([(1e12 + HEALPIX_id * 1e8 + ii).astype(
    'int').astype('str').zfill(12) for ii in n.arange(N_targets)])

t = Table()

t['NAME'] = Column(name_array, dtype=str)
t['RA'] = Column(ra_array, unit='degree', dtype=n.float64)
t['DEC'] = Column(dec_array, unit='degree', dtype=n.float64)
t['PMRA'] = Column(n.zeros(N_targets), unit='mas/yr', dtype=n.float32)
t['PMDEC'] = Column(n.zeros(N_targets), unit='mas/yr', dtype=n.float32)
t['EPOCH'] = Column(2000 * n.ones(N_targets), unit='yr', dtype=n.float32)
t['RESOLUTION'] = Column(n.ones(N_targets), unit='', dtype=n.int16)

s1_name = 'cluster_BCG'
s2_name = 'cluster_redGAL'
s3_name = 'filament_GAL'

subsurvey_name = n.zeros_like(name_array).astype('str')
subsurvey_name[:] = s1_name
subsurvey_name[subsurvey_id == 2] = s3_name
subsurvey_name[subsurvey_id == 3] = s2_name

t['SUBSURVEY'] = Column(subsurvey_name, unit='', dtype=str)

priority = 100 * n.ones(N_targets)
priority[subsurvey_id == 2] = 90
priority[subsurvey_id == 3] = 80
t['PRIORITY'] = Column(priority, unit='', dtype=n.int16)


template_names = n.zeros_like(z_array).astype('U100')

z_all = n.hstack((0., n.arange(0.3, 3., 0.2), 3.5, 4.5, 6.))
#z_all = n.hstack(( 0.0, n.arange(0.3, 3., 0.2), 6. ))
zmins = z_all[:-1]
zmaxs = z_all[1:]

ebv_1000 = (galactic_ebv_array * 1000).astype('int')
print('ebv_galaxy', n.min(ebv_1000), n.max(ebv_1000))
ebv_1_0 = (ebv_1000 > 1000)
ebv_0_5 = (ebv_1000 > 500) & (ebv_1000 <= 1000)
ebv_0_4 = (ebv_1000 > 400) & (ebv_1000 <= 500)
ebv_0_3 = (ebv_1000 > 300) & (ebv_1000 <= 400)
ebv_0_2 = (ebv_1000 > 200) & (ebv_1000 <= 300)
ebv_0_1 = (ebv_1000 > 100) & (ebv_1000 <= 200)
ebv_0_0 = (ebv_1000 <= 100)


def z_name(z0, z1): return "_zmin_" + str(int(10 * z0)).zfill(2) + \
    "_zmax_" + str(int(10 * z1)).zfill(2)


for z0, z1 in zip(zmins, zmaxs):
    zsel = (z_array >= z0) & (z_array < z1)
    template_names[(zsel)] = "4most_" + 'LRG' + \
        z_name(z0, z1) + '_EBV_0_01.fits'
    template_names[(zsel) & (ebv_0_0)] = "4most_" + 'LRG' + \
        z_name(z0, z1) + '_EBV_0_01.fits'
    template_names[(zsel) & (ebv_0_1)] = "4most_" + \
        'LRG' + z_name(z0, z1) + '_EBV_0_1.fits'
    template_names[(zsel) & (ebv_0_2)] = "4most_" + \
        'LRG' + z_name(z0, z1) + '_EBV_0_2.fits'
    template_names[(zsel) & (ebv_0_3)] = "4most_" + \
        'LRG' + z_name(z0, z1) + '_EBV_0_3.fits'
    template_names[(zsel) & (ebv_0_4)] = "4most_" + \
        'LRG' + z_name(z0, z1) + '_EBV_0_4.fits'
    template_names[(zsel) & (ebv_0_5)] = "4most_" + \
        'LRG' + z_name(z0, z1) + '_EBV_0_5.fits'
    template_names[(zsel) & (ebv_1_0)] = "4most_" + \
        'LRG' + z_name(z0, z1) + '_EBV_1_0.fits'

# sanity checks

tpls = sorted(n.array(list(set(template_names))))
print(tpls)
print('N templates used', len(tpls))
bad = (template_names == '0.0')
print(len(bad.nonzero()[0]))

N_all = len(bad)

t['TEMPLATE'] = Column(template_names, unit='', dtype=str)

ruleset_array = n.zeros_like(z_array).astype('U20')
ruleset_array[ruleset_array == "0.0"] = "RedGAL"
t['RULESET'] = Column(ruleset_array, unit='', dtype=str)

print(template_names, ruleset_array)
print(template_names.shape, ruleset_array.shape)
print(template_names.dtype, ruleset_array.dtype)


t['REDSHIFT_ESTIMATE'] = Column(z_array, unit='', dtype=n.float32)
t['REDSHIFT_ERROR'] = Column(
    n.ones(N_targets) *
    0.00001,
    unit='',
    dtype=n.float32)


t['EXTENT_PARAMETER'] = Column(n.ones(N_targets), unit='', dtype=n.float32)
t['REDDENING'] = Column(galactic_ebv_array, unit='mag', dtype=n.float32)
t['DATE_EARLIEST'] = Column(
    n.ones(N_targets) *
    59215.,
    unit='',
    dtype=n.float64)
t['DATE_LATEST'] = Column(n.ones(N_targets) * 66520, unit='', dtype=n.float64)

magnitude_4fs = get_qty('galaxy_mag_r')

# mass size relation
# https://arxiv.org/pdf/1411.6355.pdf
# Table 2 and 3. r mag line for sersic selection
# equation 2
# radius in kpc


def re_dev(M_star): return 0.16 * (M_star)**(0.1) * \
    (1 + M_star / (2.42 * 10**(10)))**(0.76 - 0.1)


def re_exp(M_star): return 0.08 * (M_star)**(0.16) * \
    (1 + M_star / (17.1 * 10**(10)))**(0.81 - 0.16)


radius_kpc = re_dev(10**stellar_mass)
radius_arcsec = cosmo.arcsec_per_kpc_proper(z_array).value * radius_kpc

# surface brightness profiles
# http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
b4 = 7.669
b1 = 1.678


def f_14_dev(r12): return gammainc(8, b4 * (0.7 / r12)**(1. / 6.))


def f_14_exp(r12): return gammainc(
    2, b1 * (0.7 / r12)**(1. / 2.))  # From Raichoor 2017


def f_14_test(r12, nn): return gammainc(2, b1 * (0.7 / r12)**(1. / nn))


frac = f_14_dev(radius_arcsec)

flux_fiber = frac * 10**((magnitude_4fs + 48.6) / -2.5)
magnitude_4fs2 = -2.5 * n.log10(flux_fiber) - 48.6

mag_type = n.zeros(N_targets).astype('U10')
mag_type[mag_type == "0.0"] = "SDSS_r_AB"

t['MAG'] = Column(magnitude_4fs2, unit='mag', dtype=n.float32)
t['MAG_ERR'] = Column(
    n.ones_like(magnitude_4fs2) *
    0.1,
    unit='',
    dtype=n.float32)
t['MAG_TYPE'] = Column(mag_type, unit='', dtype=n.str)


print(path_2_4MOST_catalog)
if os.path.isfile(path_2_4MOST_catalog):
    os.system("rm " + path_2_4MOST_catalog)

t.write(path_2_4MOST_catalog, format='fits')
