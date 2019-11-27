"""
What it does
------------

Plots the cluster model

References
----------

Command to run
--------------

python3 004_5_plot_clusters.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

"""

import glob
import sys
from astropy_healpix import healpy
import os
import time
import matplotlib.pyplot as p
import matplotlib
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


matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})

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
path_2_CLU_SAT_RS_catalog = os.path.join(test_dir, env + '_eRO_CLU_SAT_RS.fit')

fig_dir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'figures',
    env,
    'clusters',
)
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)

path_2_spiders_host_file = os.path.join(
    os.environ['HOME'],
    'data/spiders/cluster',
    'validatedclusters_catalogue_2018-12-04_version_round1-v1_Xmass1-v1.fits.gz')
spiders_host = fits.open(path_2_spiders_host_file)[1].data

path_2_spiders_file = os.path.join(
    os.environ['HOME'],
    'data/spiders/cluster',
    'validatedclusters_catalogue_2018-12-04_version_round1-v1_Xmass1-v1-flat.fits')
spiders_gal = fits.open(path_2_spiders_file)[1].data

# red sequence DATA
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

# opti0nal red sequence figure

fig_out = os.path.join(red_seq_dir, 'red_sequence.png')

p.figure(1, (6., 6.))

p.errorbar(redshift_RS, gr_RS, xerr=0.05, yerr=gr_sigma_RS, label='g-r')
p.errorbar(redshift_RS, ri_RS, xerr=0.05, yerr=ri_sigma_RS, label='r-i')
p.errorbar(redshift_RS, iz_RS, xerr=0.05, yerr=iz_sigma_RS, label='i-z')

p.legend(frameon=False)
# p.title('z='+str(z_cluster))
p.xlabel('redshift')
p.ylabel('color')
p.grid()
# p.ylim((0,1.1))
# p.yscale('log')
# p.xscale('log')
p.savefig(fig_out)
p.clf()

# opens mock catalogs
# clusters
hdu_host = fits.open(path_2_CLU_catalog)
# galaxies
hdu_clu = fits.open(path_2_CLU_SAT_RS_catalog)  # , mode='update')
x = hdu_clu[1].data['comoving_distance_to_cluster_in_rvir']
is_quiescent = hdu_clu[1].data['is_quiescent']
zr_CLU = hdu_clu[1].data['redshift_R']
richness = hdu_clu[1].data['richness']
#log_sfr = hdu_clu[1].data['galaxy_star_formation_rate']

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

# schechter functio, Eq. 8,9.


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


# TB implemented to modify this column :
mag_r = hdu_clu[1].data['galaxy_mag_r']
mag_abs_r = hdu_clu[1].data['galaxy_mag_abs_r']
DM = mag_r - mag_abs_r

vmax = hdu_clu[1].data['HALO_vmax']

mag_abs_r_HAM = n.zeros_like(mag_abs_r)

m_star_p1 = -22 + 2
# galaxies brigther than M star + 1, -21
bright_red_galaxies = (mag_abs_r <= m_star_p1) & (is_quiescent) & (x <= 1)
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

print(time.time() - t0)
z_bins = n.arange(0., n.max(zr_CLU), 0.1)
z_cluster = hdu_host[1].data['redshift_R']
for zmin, zmax in zip(z_bins, z_bins + 0.1):
    print(zmin, '<z<', zmax)
    # define host clusters, area in a redshift bin
    sel_clusters = (z_cluster >= zmin) & (z_cluster < zmax)
    z_cluster_mean = n.mean(z_cluster[sel_clusters])
    print('z_cluster_mean', z_cluster_mean)
    area_Mpc2 = n.sum(
        n.pi * (hdu_host[1].data['HALO_rvir'][sel_clusters] / 1000)**2)
    print('area_Mpc2', area_Mpc2)
    # properties of galaxy members
    sel = (zr_CLU >= zmin) & (zr_CLU < zmax)
    richness_in_cluster = n.mean(richness[sel])
    #magnitudes = mag_abs_r[sel]
    N_in_cluster = len(richness[sel])
    rds = norm.rvs(loc=0, scale=0.15, size=N_in_cluster)
    print('N_in_cluster', N_in_cluster)
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
    inter = mag_abs_r_HAM[sel]
    print('inter', inter, inter)
    inter[id_sort_vmax] = mags_2_assign
    mag_abs_r_HAM[sel] = inter
    #logPhi_evol(zr_CLU, richness)
    print(time.time() - t0)
    # now the SPIDERS DATA
    sel_spiders_host = (
        spiders_host['NMEMBERS'] >= 6) & (
        spiders_host['CLUZSPEC'] >= zmin) & (
            spiders_host['CLUZSPEC'] < zmax)
    z_cluster_spiders = spiders_host['CLUZSPEC'][sel_spiders_host]
    z_cluster_spiders_mean = n.mean(spiders_host['CLUZSPEC'][sel_spiders_host])
    conversion = (
        1 /
        cosmo.arcsec_per_kpc_proper(z_cluster_spiders).to(
            u.deg /
            u.Mpc)).value
    radius_Mpc_spiders = spiders_host['R200C_DEG'][sel_spiders_host] * conversion
    area_Mpc2_spiders = n.sum(n.pi * (radius_Mpc_spiders)**2)
    sel_spiders = (
        spiders_gal['ISMEMBER'] == 1) & (
        spiders_gal['CLUZSPEC'] >= zmin) & (
            spiders_gal['CLUZSPEC'] < zmax)
    abs_mag_i_spiders = spiders_gal['CMODELMAG_I'][sel_spiders] - \
        cosmo.distmod(spiders_gal['IDLSPEC1D_Z'][sel_spiders]).value
    # figure of LF
    fig_out = os.path.join(fig_dir, 'cluster_galaxy_LF_z_' +
                           str(n.round(z_cluster_mean, 3)) + '.png')
    mdex = 0.25
    mbins = n.arange(-30, -10, mdex)
    p.figure(1, (6., 5.5))
    p.axes([0.15, 0.15, 0.8, 0.8])
    # p.tight_layout()
    #p.hist(mag_abs_r[sel]    , lw=3, weights=n.ones(N_in_cluster)/mdex/area_Mpc2, bins=mbins, histtype='step', rasterized=True, label='magneticum, r')
    out = p.hist(
        mag_abs_r_HAM[sel],
        lw=3,
        weights=n.ones(N_in_cluster) /
        mdex /
        area_Mpc2,
        bins=mbins,
        histtype='step',
        rasterized=True,
        label='Mock R18, r')
    p.hist(mag_abs_r_HAM[sel & (is_quiescent == False) & (x <= 1)],
           lw=3,
           weights=n.ones_like(mag_abs_r_HAM[sel & (is_quiescent == False) & (x <= 1)]) / mdex / area_Mpc2,
           bins=mbins,
           histtype='step',
           rasterized=True,
           label='Mock, star-forming')
    p.hist(mag_abs_r_HAM[sel & is_quiescent & (x <= 1)],
           lw=3,
           weights=n.ones_like(mag_abs_r_HAM[sel & is_quiescent & (x <= 1)]) / mdex / area_Mpc2,
           bins=mbins,
           histtype='step',
           rasterized=True,
           label='Mock quiescent')
    p.hist(mag_abs_r_HAM[sel & is_quiescent & (x <= .1)],
           lw=3,
           weights=n.ones_like(mag_abs_r_HAM[sel & is_quiescent & (x <= 0.1)]) / mdex / area_Mpc2,
           bins=mbins,
           histtype='step',
           rasterized=True,
           label='Mock quiescent r<0.1rvir',
           ls='dashed')
    p.hist(
        abs_mag_i_spiders,
        lw=2,
        weights=n.ones(
            len(abs_mag_i_spiders)) /
        mdex /
        area_Mpc2_spiders,
        bins=mbins,
        histtype='step',
        rasterized=True,
        label=r'SPIDERS, i, $N_{mem}>6$')
    p.title(r'$\bar{z}$=' + str(n.round(z_cluster_mean, 3)))
    p.ylabel(r'$\Phi=N_{gal} mag^{-1} Mpc^{-2}$')
    p.xlabel('Absolute magnitude')
    p.yscale('log')
    p.xlim((n.min(out[1][1:][out[0] > 0]) - mdex,
            n.max(out[1][1:][out[0] > 0]) + mdex))
    p.ylim((1e-3, 20))
    p.legend(loc=4, frameon=False)
    p.grid()
    p.savefig(fig_out)
    p.clf()

fig_out = os.path.join(fig_dir, 'cluster_galaxy_LF-mag-mag.png')

p.figure(1, (6., 5.5))
p.tight_layout()
p.plot(mag_abs_r, mag_abs_r_HAM, 'k,', rasterized=True)
p.plot(mrange, mrange, 'm--', lw=2)
# p.title(baseName)
# p.yscale('log')
p.xlabel('mag MAGNETICUM')
p.ylabel('mag HAM')
p.grid()
p.savefig(fig_out)
p.clf()

fig_out = os.path.join(fig_dir, 'cluster_galaxy_z-mag-comparison.png')

# &(spiders_gal['CLUZSPEC']>=zmin)&(spiders_gal['CLUZSPEC']<zmax)
sel_spiders = (spiders_gal['ISMEMBER'] == 1)
abs_mag_i_spiders = spiders_gal['CMODELMAG_I'][sel_spiders] - \
    cosmo.distmod(spiders_gal['IDLSPEC1D_Z'][sel_spiders]).value
z_spiders = spiders_gal['IDLSPEC1D_Z'][sel_spiders]

p.figure(1, (6., 5.5))
p.tight_layout()
p.plot(zr_CLU, mag_abs_r, 'r,', rasterized=True, label='magneticum')
p.plot(zr_CLU, mag_abs_r_HAM, 'k,', rasterized=True, label='HAM Ricci 18')
p.plot(z_spiders, abs_mag_i_spiders, 'b,', rasterized=True, label='SPIDERS')
p.legend(loc=0)
# p.title(baseName)
# p.yscale('log')
p.xlabel('redshift')
p.ylabel('mag')
p.grid()
p.savefig(fig_out)
p.clf()


fig_out = os.path.join(fig_dir, 'cluster_galaxy_LF.png')

p.figure(1, (6., 5.5))
p.tight_layout()
p.plot(zr_CLU, mag_abs_r, 'r,', rasterized=True, label='magneticum')
p.plot(zr_CLU, mag_abs_r_HAM, 'k,', rasterized=True, label='HAM Ricci 18')
# p.title(baseName)
# p.yscale('log')
p.xlabel('redshift')
p.ylabel('mag')
p.grid()
p.savefig(fig_out)
p.clf()

fig_out = os.path.join(fig_dir, 'cluster_vmax-magnitude.png')

p.figure(1, (6., 5.5))
p.tight_layout()
p.plot(
    hdu_clu[1].data['HALO_vmax'],
    mag_abs_r,
    'r,',
    rasterized=True,
    label='magneticum')
p.plot(hdu_clu[1].data['HALO_vmax'], mag_abs_r_HAM,
       'k,', rasterized=True, label='HAM Ricci 18')
# p.title(baseName)
# p.yscale('log')
p.xlabel('V max')
p.ylabel('mag')
p.xscale('log')
p.grid()
p.savefig(fig_out)
p.clf()

test2 = n.unique(
    hdu_clu[1].data['HALO_host_id'],
    return_counts=True,
    return_index=True,
    return_inverse=True)
print('time elapsed', time.time() - t0)
print('N clusers', len(test2[0]))
for i_clu in test2[0][202337:]:
    sel = (test2[2] == i_clu)
    z_in_cluster = zr_CLU[sel]
    magnitudes = mag_abs_r[sel]
    vmax_in_cluster = vmax[sel]
    richness_in_cluster = richness[sel][0]
    N_in_cluster = len(z_in_cluster)
    area_Mpc2 = n.pi * (hdu_host[1].data['HALO_rvir'][i_clu] / 1000)**2
    z_cluster = hdu_host[1].data['redshift_R'][i_clu]
    mrange = n.arange(-30, -15, 0.01)
    total_number = n.cumsum(
        Schechter_M_z(
            mrange,
            z_cluster,
            richness_in_cluster) *
        area_Mpc2 *
        0.01)
    x_itp = n.hstack(
        (total_number[total_number > 0.1], total_number[total_number > 0.1][-1] + 1000))
    y_itp = n.hstack((mrange[total_number > 0.1],
                      mrange[total_number > 0.1][-1]))
    itp_clu = interp1d(x_itp, y_itp)
    mags_2_assign = itp_clu(n.arange(N_in_cluster) + 1)
    id_sort_vmax = n.argsort(vmax_in_cluster)
    # id=0 the lowest vmax
    inter = mag_abs_r_HAM[sel]
    inter[id_sort_vmax] = mags_2_assign
    mag_abs_r_HAM[sel] = inter
    #logPhi_evol(zr_CLU, richness)
print('time elapsed', time.time() - t0)

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

# create colors for the galaxies

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

# Add uncertainties on the magnitudes and colors
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

NSIDE = 2048
# ordering='RING'

#ra,dec = healpy.pix2ang(NSIDE, n.arange(healpy.nside2npix(NSIDE)),nest=False, lonlat=True)

sdss_g_depth = healpy.fitsfunc.read_map(depth_sdss_g)
sdss_r_depth = healpy.fitsfunc.read_map(depth_sdss_r)
sdss_i_depth = healpy.fitsfunc.read_map(depth_sdss_i)
sdss_z_depth = healpy.fitsfunc.read_map(depth_sdss_z)

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

sdss_g_temp = mag_abs_r_HAM + DM2 + gr_all
sdss_r_temp = mag_abs_r_HAM + DM2
sdss_i_temp = mag_abs_r_HAM + DM2 - ri_all
sdss_z_temp = mag_abs_r_HAM + DM2 - ri_all - iz_all

#
sdss_g_err = 2.5 * 10**(-0.4 * (depth_10_sigma_g - sdss_g_temp)) / n.log(10)
sdss_r_err = 2.5 * 10**(-0.4 * (depth_10_sigma_r - sdss_r_temp)) / n.log(10)
sdss_i_err = 2.5 * 10**(-0.4 * (depth_10_sigma_i - sdss_i_temp)) / n.log(10)
sdss_z_err = 2.5 * 10**(-0.4 * (depth_10_sigma_z - sdss_z_temp)) / n.log(10)

sdss_g_err[sdss_g_err == n.inf] = 0
sdss_r_err[sdss_r_err == n.inf] = 0
sdss_i_err[sdss_i_err == n.inf] = 0
sdss_z_err[sdss_z_err == n.inf] = 0

N_gal = len(sdss_g_temp)
# random variable in a norm distribution with sigma^2 = 1/40
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
