"""
What it does
------------

Computes magnitudes for galaxies in clusters following the luminosity function of Ricci et al. 2018

Adds color to red sequence galaxies
Model taken from redmapper trained on SDSS-IV/SPIDERS by J. Ider Chitham

Make a realization of galaxies following SDSS depth maps

outputs to a single file: $ENV_eRO_CLU_SAT_RS.fit

References
----------

 * Ricci et al. 2018
 * Rykoff et al. 2013
 * Clerc et al. 2016

Command to run
--------------

python3 004_4_red_sequence.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

red sequence models

luminosity function

"""

import glob
import sys
from astropy_healpix import healpy
import os
import time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import healpy as hp
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


env = sys.argv[1]

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
test_dir = os.path.join(os.environ[env], 'hlists', 'fits')
path_2_CLU_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')
path_2_CLU_SAT_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT.fit')
# output catalog
path_2_out_file = os.path.join(test_dir, env + '_eRO_CLU_SAT_RS.fit')

fig_dir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'figures',
    env,
    'clusters',
)
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)

# ==============================
# READ red sequence DATA
# ==============================
red_seq_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'redSequence')
path_2_f1 = os.path.join(red_seq_dir, 'sdss_cal_iter1_pars.fit')
path_2_f2 = os.path.join(red_seq_dir, 'sdss_cal_zspec_redgals_model.fit')

f1 = fits.open(path_2_f1)
f2 = fits.open(path_2_f2)

redshift_RS = f2[1].data['nodes'][0]
ug_RS = f2[1].data['meancol'].T[0].T[0]
gr_RS = f2[1].data['meancol'].T[1].T[0]
ri_RS = f2[1].data['meancol'].T[2].T[0]
iz_RS = f2[1].data['meancol'].T[3].T[0]

ug_sigma_RS = f2[1].data['meancol_scatter'].T[0].T[0]
gr_sigma_RS = f2[1].data['meancol_scatter'].T[1].T[0]
ri_sigma_RS = f2[1].data['meancol_scatter'].T[2].T[0]
iz_sigma_RS = f2[1].data['meancol_scatter'].T[3].T[0]

ug_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
                     n.hstack((ug_RS[0], ug_RS, ug_RS[-1])))
gr_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
                     n.hstack((gr_RS[0], gr_RS, gr_RS[-1])))
ri_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
                     n.hstack((ri_RS[0], ri_RS, ri_RS[-1])))
iz_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)),
                     n.hstack((iz_RS[0], iz_RS, iz_RS[-1])))

ug_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
    (ug_sigma_RS[0], ug_sigma_RS, ug_sigma_RS[-1])))
gr_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
    (gr_sigma_RS[0], gr_sigma_RS, gr_sigma_RS[-1])))
ri_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
    (ri_sigma_RS[0], ri_sigma_RS, ri_sigma_RS[-1])))
iz_sigma_RS_itp = interp1d(n.hstack((0., redshift_RS, 2.0)), n.hstack(
    (iz_sigma_RS[0], iz_sigma_RS, iz_sigma_RS[-1])))

# ==============================
# Luminosity function model
# ==============================
# abundance matching for magnitudes based on Ricci et al. 2018. Tables
# 3,5, Eq. 8,9,10
logPhi_z, alpha_z, M_s_z = n.loadtxt(os.path.join(
    os.environ['GIT_AGN_MOCK'], 'data', 'cluster_galaxies', 'ricci_2018_table3.ascii'))

# redshift evolution function, Eq. 10


def fun(redshift, richness, a, b, c): return a * \
    n.log10(1 + redshift) + b * n.log10(richness) + c


def logPhi_evol(
    redshift,
    richness): return fun(
        redshift,
        richness,
        logPhi_z[0],
        logPhi_z[1],
    logPhi_z[2])


def alpha_evol(
    redshift,
    richness): return fun(
        redshift,
        richness,
        alpha_z[0],
        alpha_z[1],
    alpha_z[2])


def M_s_evol(
    redshift,
    richness): return fun(
        redshift,
        richness,
        M_s_z[0],
        M_s_z[1],
    M_s_z[2])

# schechter function, Eq. 8,9.


def Schechter_L(L, phi_s, L_s, alpha): return phi_s * \
    (L / L_s)**alpha * n.e**(-L / L_s) / L_s


def Schechter_M(M, phi_s, M_s, alpha): return 0.4 * n.log(10.) * phi_s * \
    10**(0.4 * (M_s - M) * (alpha + 1)) * n.e**(-10**(0.4 * (M_s - M)))


def Schechter_M_z(M,
                  redshift,
                  richness): return 0.4 * n.log(10.) * 10**logPhi_evol(redshift,
                                                                       richness) * 10**(0.4 * (M_s_evol(redshift,
                                                                                                        richness) - M) * (alpha_evol(redshift,
                                                                                                                                     richness) + 1)) * n.e**(-10**(0.4 * (M_s_evol(redshift,
                                                                                                                                                                                   richness) - M)))


# ==============================
# READ mock catalogs
# ==============================
# clusters
hdu_host = fits.open(path_2_CLU_catalog)

# galaxies
hdu_clu = fits.open(path_2_CLU_SAT_catalog)  # , mode='update')
x = hdu_clu[1].data['comoving_distance_to_cluster_in_rvir']
is_quiescent = hdu_clu[1].data['is_quiescent']
zr_CLU = hdu_clu[1].data['redshift_R']
vmax = hdu_clu[1].data['HALO_vmax']

N_gal = len(zr_CLU)
# initializes the array magnitude to be computed to 0s
mag_abs_r_HAM_ini = hdu_clu[1].data['galaxy_mag_abs_r']
#mag_abs_r_HAM = hdu_clu[1].data['galaxy_mag_abs_r']
mag_abs_r_HAM = n.zeros_like(vmax)

# ==============================
## COMPUTES THE RICHNESS a first time ##
# ==============================
# select galaxies brigther than M star + 2, -20
m_star_p1 = -22 + 2
bright_red_galaxies = (
    mag_abs_r_HAM_ini <= m_star_p1) & (is_quiescent) & (
        x <= 1)
faint_galaxies = (~bright_red_galaxies)
# richness estimator, for red galaxies around the cluster such as m <
# m_star + 2
test3 = n.unique(
    hdu_clu[1].data['HALO_host_id'][bright_red_galaxies],
    return_counts=True,
    return_index=True,
    return_inverse=True)
#N_occurences_3 = test3[3][test3[2]]
richness_per_cluster = interp1d(test3[0], test3[3])
richness = richness_per_cluster(hdu_clu[1].data['HALO_host_id'])
richness_col = fits.Column(
    name='richness',
    format='K',
    unit='',
    array=richness)
richness1 = n.copy(richness)
print('Richness 1st estimate (M star + 2): ', richness[:10], time.time() - t0)

print('=============================================')
print('=============================================')
print('Abundance matching')
# ==============================
## ABUNDANCE MATCHING ##
# ==============================
# assigns magnitudes to the galaxies in the clusters
z_bins = n.arange(0., n.max(zr_CLU), 0.1)
z_cluster = hdu_host[1].data['redshift_R']
for zmin, zmax in zip(z_bins, z_bins + 0.1):
    print(zmin, '<z<', zmax)
    # select host clusters, in a redshift bin, compute area in Mpc2
    sel_clusters = (z_cluster >= zmin) & (z_cluster < zmax)
    z_cluster_mean = n.mean(z_cluster[sel_clusters])
    print('z_cluster_mean', z_cluster_mean,
          'N clusters=', len(z_cluster[sel_clusters]))
    area_Mpc2 = n.sum(
        n.pi * (hdu_host[1].data['HALO_rvir'][sel_clusters] / 1000)**2)
    print('area_Mpc2', area_Mpc2)
    # select member galaxies
    sel = (zr_CLU >= zmin) & (zr_CLU < zmax)
    richness_in_cluster = n.mean(richness[sel])
    N_in_cluster = len(richness[sel])
    print('N galaxies _in_cluster', N_in_cluster)
    # scatter in the abundance matching relation: random Gaussian number,
    # scale=0.15, loc=0
    rds = norm.rvs(loc=0, scale=0.15, size=N_in_cluster)
    vmax_in_cluster = 10**(n.log10(vmax[sel]) + rds)
    mrange = n.arange(-30, -10, 0.01)
    total_number = n.cumsum(
        Schechter_M_z(
            mrange,
            z_cluster_mean,
            richness_in_cluster) *
        area_Mpc2 *
        0.01)
    print('total_number', total_number)
    x_itp = n.hstack(
        (total_number[total_number > 0.1], total_number[total_number > 0.1][-1] * 10))
    y_itp = n.hstack((mrange[total_number > 0.1],
                      mrange[total_number > 0.1][-1]))
    itp_clu = interp1d(x_itp, y_itp)
    mags_2_assign = itp_clu(n.arange(N_in_cluster) + 1)
    print('mags_2_assign', mags_2_assign, mags_2_assign.shape)
    id_sort_vmax = n.argsort(vmax_in_cluster)[::-1]
    # id=0 the lowest vmax
    inter = n.zeros_like(mag_abs_r_HAM[sel])
    print('inter', inter)
    inter[id_sort_vmax] = mags_2_assign
    print('inter', inter)
    mag_abs_r_HAM[sel] = inter
    print('mag_abs_r_HAM[sel]', mag_abs_r_HAM[sel])
    #logPhi_evol(zr_CLU, richness)
    print('dt=', time.time() - t0)


z_bins_2 = n.arange(0., n.max(zr_CLU) + 0.1, 0.001)
DM2_itp = interp1d(z_bins_2, cosmo.distmod(z_bins_2).value)
DM2 = DM2_itp(zr_CLU)
rmag_abs_ham_col = fits.Column(
    name='galaxy_mag_abs_r',
    format='D',
    unit='mag',
    array=mag_abs_r_HAM)
rmag_ham_col = fits.Column(
    name='galaxy_mag_r',
    format='D',
    unit='mag',
    array=mag_abs_r_HAM + DM2)
print(
    'mean(initial magnitude - new magnitude)',
    n.mean(
        mag_abs_r_HAM_ini -
        mag_abs_r_HAM),
    time.time() -
    t0)
print('HAM r magnitudes, time elapsed', mag_abs_r_HAM[:10], time.time() - t0)

# ==============================
## COMPUTES THE RICHNESS a second time ##
# ==============================
bright_red_galaxies = (mag_abs_r_HAM <= m_star_p1) & (is_quiescent) & (x <= 1)
faint_galaxies = (~bright_red_galaxies)
test3 = n.unique(
    hdu_clu[1].data['HALO_host_id'][bright_red_galaxies],
    return_counts=True,
    return_index=True,
    return_inverse=True)
richness_per_cluster = interp1d(test3[0], test3[3])
richness = richness_per_cluster(hdu_clu[1].data['HALO_host_id'])
richness_col = fits.Column(
    name='richness',
    format='K',
    unit='',
    array=richness)
print('Richness 2nd estimate (M star + 2): ', richness[:10], time.time() - t0)
print('mean richness difference', n.mean(richness1 - richness))

# ==============================
## COMPUTES THE COLOR ##
# ==============================
rds = norm.rvs(loc=0, scale=0.001, size=len(zr_CLU[is_quiescent]))

# + gr_sigma_RS_itp (zr_CLU[is_quiescent]) * rds / 10.
gr = gr_RS_itp(zr_CLU[is_quiescent]) + rds
# + ri_sigma_RS_itp (zr_CLU[is_quiescent]) * rds / 10.
ri = ri_RS_itp(zr_CLU[is_quiescent]) + rds
# + iz_sigma_RS_itp (zr_CLU[is_quiescent]) * rds / 10.
iz = iz_RS_itp(zr_CLU[is_quiescent]) + rds

gr_all = n.zeros_like(zr_CLU)
ri_all = n.zeros_like(zr_CLU)
iz_all = n.zeros_like(zr_CLU)
gr_all[is_quiescent] = gr
ri_all[is_quiescent] = ri
iz_all[is_quiescent] = iz

gr_col = fits.Column(name='galaxy_gr', format='D', unit='mag', array=gr_all)
ri_col = fits.Column(name='galaxy_ri', format='D', unit='mag', array=ri_all)
iz_col = fits.Column(name='galaxy_iz', format='D', unit='mag', array=iz_all)
print('colors assigned',
      gr_all[:10],
      ri_all[:10],
      iz_all[:10],
      time.time() - t0)

# ==============================
## COMPUTES SDSS magnitudes ##
# ==============================
# Add uncertainties on the magnitudes and colors following SDSS depth maps
# reads the maps
sdss_depth_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'data', 'sdss_depth')
depth_sdss_g = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_g_model_10sigma.fits.gz')
depth_sdss_i = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_i_model_10sigma.fits.gz')
depth_sdss_r = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_r_model_10sigma.fits.gz')
depth_sdss_z = os.path.join(
    sdss_depth_dir,
    'sdss_dr8_nodered_nside2048_z_model_10sigma.fits.gz')
#
NSIDE = 2048
sdss_g_depth = hp.fitsfunc.read_map(depth_sdss_g)
sdss_r_depth = hp.fitsfunc.read_map(depth_sdss_r)
sdss_i_depth = hp.fitsfunc.read_map(depth_sdss_i)
sdss_z_depth = hp.fitsfunc.read_map(depth_sdss_z)
# retrieves the depths
pixels = healpy.ang2pix(
    NSIDE,
    hdu_clu[1].data['ra'],
    hdu_clu[1].data['dec'],
    nest=False,
    lonlat=True)
depth_10_sigma_g = sdss_g_depth[pixels]
depth_10_sigma_r = sdss_r_depth[pixels]
depth_10_sigma_i = sdss_i_depth[pixels]
depth_10_sigma_z = sdss_z_depth[pixels]
# computes the magnitudes
sdss_g_temp = mag_abs_r_HAM + DM2 + gr_all
sdss_r_temp = mag_abs_r_HAM + DM2
sdss_i_temp = mag_abs_r_HAM + DM2 - ri_all
sdss_z_temp = mag_abs_r_HAM + DM2 - ri_all - iz_all
# their uncertainty
sdss_g_err = 2.5 * 10**(-0.4 * (depth_10_sigma_g - sdss_g_temp)) / n.log(10)
sdss_r_err = 2.5 * 10**(-0.4 * (depth_10_sigma_r - sdss_r_temp)) / n.log(10)
sdss_i_err = 2.5 * 10**(-0.4 * (depth_10_sigma_i - sdss_i_temp)) / n.log(10)
sdss_z_err = 2.5 * 10**(-0.4 * (depth_10_sigma_z - sdss_z_temp)) / n.log(10)
# assign uncertainty=0 outside of the SDSS footprint
sdss_g_err[sdss_g_err == n.inf] = 0
sdss_r_err[sdss_r_err == n.inf] = 0
sdss_i_err[sdss_i_err == n.inf] = 0
sdss_z_err[sdss_z_err == n.inf] = 0
print('SDSS magnitude errors computed',
      sdss_g_err[:10],
      sdss_r_err[:10],
      sdss_i_err[:10],
      sdss_z_err[:10],
      time.time() - t0)

# random variable in a norm distribution with sigma = 1

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_g_col = fits.Column(
    name='sdss_g',
    format='D',
    unit='mag',
    array=sdss_g_temp +
    rds *
    sdss_g_err)

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_r_col = fits.Column(
    name='sdss_r',
    format='D',
    unit='mag',
    array=sdss_r_temp +
    rds *
    sdss_r_err)

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_i_col = fits.Column(
    name='sdss_i',
    format='D',
    unit='mag',
    array=sdss_i_temp +
    rds *
    sdss_i_err)

rds = norm.rvs(loc=0, scale=1, size=N_gal)
sdss_z_col = fits.Column(
    name='sdss_z',
    format='D',
    unit='mag',
    array=sdss_z_temp +
    rds *
    sdss_z_err)


sdss_g_err_col = fits.Column(
    name='sdss_g_err',
    format='D',
    unit='mag',
    array=sdss_g_err)
sdss_r_err_col = fits.Column(
    name='sdss_r_err',
    format='D',
    unit='mag',
    array=sdss_r_err)
sdss_i_err_col = fits.Column(
    name='sdss_i_err',
    format='D',
    unit='mag',
    array=sdss_i_err)
sdss_z_err_col = fits.Column(
    name='sdss_z_err',
    format='D',
    unit='mag',
    array=sdss_z_err)


print('creates hdu, time elapsed', time.time() - t0)

all_cols = []
for name, uni, forma in zip(
        hdu_clu[1].data.columns.names, hdu_clu[1].data.columns.units, hdu_clu[1].data.columns.formats):
    if name == "galaxy_mag_abs_r" or name == "galaxy_mag_r":
        print(name)
    else:
        all_cols.append(
            fits.Column(
                name=name,
                format=forma,
                unit=uni,
                array=hdu_clu[1].data[name]))

all_cols.append(gr_col)
all_cols.append(ri_col)
all_cols.append(iz_col)
all_cols.append(richness_col)
all_cols.append(rmag_abs_ham_col)
all_cols.append(rmag_ham_col)

all_cols.append(sdss_g_col)
all_cols.append(sdss_r_col)
all_cols.append(sdss_i_col)
all_cols.append(sdss_z_col)

all_cols.append(sdss_g_err_col)
all_cols.append(sdss_r_err_col)
all_cols.append(sdss_i_err_col)
all_cols.append(sdss_z_err_col)


new_cols = fits.ColDefs(all_cols)
tbhdu = fits.BinTableHDU.from_columns(new_cols)

prihdr = fits.Header()
prihdr['author'] = 'JC'
prihdu = fits.PrimaryHDU(header=prihdr)

hdu = fits.HDUList([prihdu, tbhdu])

print('writing', time.time() - t0)

if os.path.isfile(path_2_out_file):
    os.remove(path_2_out_file)

hdu.writeto(path_2_out_file)


#test = fits.open(path_2_out_file)
# print(test[1].data.columns)
