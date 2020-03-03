"""
Adds booleans corresponding to the 4MOST COSMO tracers

Uses HOD models.

output:
f['/4MOST_COSMO_tracers/is_ELG']
f['/4MOST_COSMO_tracers/is_QSO']
f['/4MOST_COSMO_tracers/is_BG_lz']
f['/4MOST_COSMO_tracers/is_BG_hz']

example:
python3.4 lc_add_4MOST_COSMO.py MD10 all_1.00000
"""
from astropy_healpix import healpy
import astropy.constants as cc
import astropy.io.fits as fits
from scipy.special import erf
from scipy.interpolate import interp2d
import h5py    # HDF5 support
import random
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.stats import norm
import sys
import numpy as n
import os
import glob
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import time
t0 = time.time()

print('4MOST COSMO TRACERS')
print('------------------------------------------------')
print('------------------------------------------------')
cosmoMD = FlatLambdaCDM(
    H0=67.77 * u.km / u.s / u.Mpc,
    Om0=0.307115)  # , Ob0=0.048206)
h = 0.6777


# specific functions


env = 'MD04'  # sys.argv[1] # 'MD04'
baseName = 'all_1.00000'  # sys.argv[2] # "sat_0.62840"
print(env, baseName)
area = 129600. / n.pi

if env == 'MD04':
    git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'DM_SMDPL')
    test_dir = os.path.join(os.environ[env], 'fits')
    L_box = 400.0 / h

if env == 'MD10':
    git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'DM_MDPL2')
    test_dir = os.path.join(os.environ[env], 'hlist', 'fits')
    L_box = 1000.0 / h

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.fits')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
path_2_agn_file = os.path.join(test_dir, baseName + '_agn.fits')
path_2_eRO_catalog = os.path.join(test_dir, baseName + '_eRO_AGN.fit')


# redshift and
N_snap, Z_snap, A_snap, DC_max, DC_min = n.loadtxt(os.path.join(
    os.environ['GIT_AGN_MOCK'], 'python', 'DM_SMDPL', 'snap_list_with_border.txt'), unpack=True)


def get_a(path_2_in):
    alp = os.path.basename(path_2_in).split('_')[1][:-5]
    print('a=', alp)
    return float(alp)


a_snap = get_a(path_2_light_cone)
sel = (A_snap == a_snap)

z_array = n.arange(0, 20., 0.0005)
dcs = cosmoMD.comoving_distance(z_array)
dc_to_z = interp1d(dcs, z_array)
z_to_dc = interp1d(z_array, dcs)

z_start = dc_to_z(DC_min[sel])[0]
z_reach = dc_to_z(DC_max[sel])[0]
z_mean = 0.5 * (z_start + z_reach)
print(z_start, '<z<', z_reach)
vol = (cosmoMD.comoving_volume(z_reach).value -
       cosmoMD.comoving_volume(z_start).value)
DL_mean_z = (cosmoMD.luminosity_distance(z_mean).to(u.cm)).value
print('volume', vol, 'Mpc3')


# HALO File
f1 = fits.open(path_2_light_cone)
if f1[1].columns.names[0] == 'id':
    Mvir = f1[1].data['Mvir'] / h
    vmax = f1[1].data['vmax']
    M200c = n.log10(f1[1].data['M200c'] / h)

else:
    hh = 'id Mvir Rvir rs scale_of_last_MM vmax x y z vx vy vz M200c M500c b_to_a_500c c_to_a_500c Acc_Rate_1_Tdyn'
    header = hh.split()
    dc_col = {}
    for ii in range(len(header)):
        dc_col[header[ii]] = 'col' + str(ii + 1)

    Mvir = f1[1].data[dc_col['Mvir']] / h
    M200c = n.log10(f1[1].data[dc_col['M200c']] / h)
    vmax = f1[1].data[dc_col['vmax']]
f1.close()
N_halos = len(vmax)

# ordering with scatter
rds = norm.rvs(loc=0, scale=0.3, size=N_halos)
vmax_sort_id = n.argsort(vmax * 10**rds)

# COORDINATE File
f2 = h5py.File(path_2_coordinate_file, 'r')
zz = f2[1].data['redshift_R'].value
dL_cm = f2[1].data['dL'].value
galactic_NH = f2['nH'].value
galactic_ebv = f2[1].data['ebv'].value
g_lat = f2[1].data['g_lat'].value
g_lon = f2[1].data['g_lon'].value
N_galaxies = len(zz)
f2.close()

# GALAXY file
f3 = h5py.File(path_2_galaxy_file, 'r')
log_stellar_mass = f3['/galaxy/SMHMR_mass'].value  # log of the stellar mass
mag_abs_r = f3['/galaxy/mag_abs_r'].value
mag_r = f3['/galaxy/mag_r'].value
f3.close()

# AGN file
f4 = h5py.File(path_2_agn_file, 'r')
FX_hard = f4["AGN/FX_hard"].value
FX_soft = f4["AGN/FX_soft"].value
FX_soft_attenuated = f4["AGN/FX_soft_attenuated"].value
LX_hard = f4["AGN/LX_hard"].value
LX_soft = f4["AGN/LX_soft"].value
SDSS_r_AB = f4["AGN/SDSS_r_AB"].value
SDSS_r_AB_attenuated = f4["AGN/SDSS_r_AB_attenuated"].value
agn_type = f4["AGN/agn_type"].value
ids_active = f4["AGN/ids_active"].value
logNH = f4["AGN/logNH"].value
random = f4["AGN/random"].value
f4.close()
mag_r_agn = SDSS_r_AB_attenuated
print('agn file opened', time.time() - t0)

# type1 AGN
ids_type1 = ids_active[((agn_type == 11) | (
    agn_type == 12)) & (SDSS_r_AB < 23.5)]
n_agn = len(ids_active)
n_agn_t1 = len(ids_type1)

# REHASH TO SHELL BOUNDARIES


def rehash(data):
    zmin, zmax, N_p_deg2 = data
    dndz = N_p_deg2 / (zmax - zmin)
    itp = interp1d(n.hstack((0., 0.5 * (zmin + zmax))), n.hstack((0., dndz)))
    total_per_deg2_in_shell = quad(itp, z_start, z_reach)[0]
    return total_per_deg2_in_shell


N_lrg1_pdeg2 = rehash(
    n.loadtxt(
        os.path.join(
            os.environ['GIT_AGN_MOCK'],
            "data/cosmo-4most/bg.nz"),
        unpack=True))
N_lrg2_pdeg2 = rehash(
    n.loadtxt(
        os.path.join(
            os.environ['GIT_AGN_MOCK'],
            "data/cosmo-4most/lrg.nz"),
        unpack=True))
N_elg_pdeg2 = rehash(
    n.loadtxt(
        os.path.join(
            os.environ['GIT_AGN_MOCK'],
            "data/cosmo-4most/elg.nz"),
        unpack=True))
N_qso_pdeg2 = rehash(
    n.loadtxt(
        os.path.join(
            os.environ['GIT_AGN_MOCK'],
            "data/cosmo-4most/qso.nz"),
        unpack=True))

print(N_lrg1_pdeg2, N_lrg2_pdeg2, N_elg_pdeg2, N_qso_pdeg2)
lrg1_selection = (n.ones(N_halos) == 0)
lrg2_selection = (n.ones(N_halos) == 0)
qso_selection = (n.ones(N_halos) == 0)
elg_selection = (n.ones(N_halos) == 0)


# LRG1 CASE

print("LRG 1", time.time() - t0)
N_lrg1 = int(area * N_lrg1_pdeg2)
print('N_lrg1', N_lrg1, time.time() - t0)
if N_lrg1 >= 1:
    ids_lrg1 = vmax_sort_id[-N_lrg1:]
    print('selected', len(ids_lrg1), 'wanted', N_lrg1, 'out of ', N_halos)


print("LRG 2", time.time() - t0)
N_lrg2 = int(area * N_lrg2_pdeg2)
print('N_lrg2', N_lrg2, time.time() - t0)
if N_lrg2 >= 1:
    # LRG1: select on Ms
    ids_lrg2 = vmax_sort_id[-N_lrg2 - N_lrg1:-N_lrg1]
    print('selected', len(ids_lrg2), 'wanted', N_lrg2, 'out of ', N_halos)


# ELG parameters
# ELG select on Mvir
mh_mean = 11.0 + z_mean
mh_scatter = 0.3

print("==============")
print("ELG", time.time() - t0)
N_elg = int(area * N_elg_pdeg2)
# if N_elg > 1 :
rds = norm.rvs(loc=0, scale=mh_scatter, size=N_halos)
M200c_scattered = M200c + rds - mh_mean
select_init = (
    M200c_scattered > -
    mh_scatter *
    1) & (
        M200c_scattered < mh_scatter *
    1)
N_avail = len(select_init.nonzero()[0])
if N_avail > N_elg:
    ids_elg = n.arange(N_halos)[select_init]
    RD = n.random.random(len(ids_elg))
    sel = (RD < N_elg * 1. / N_avail)
    ids_elg = ids_elg[sel]
else:
    select_init = (
        M200c_scattered > -
        mh_scatter *
        2) & (
        M200c_scattered < mh_scatter *
        2)
    ids_elg = n.arange(N_halos)[select_init]
    RD = n.random.random(len(ids_elg))
    sel = (RD < N_elg * 1. / N_avail)
    ids_elg = ids_elg[sel]

print("==============")
# print("==============")
print("QSO", time.time() - t0)
N_qso = int(area * N_qso_pdeg2)  # 2
if N_qso > 1 and N_qso < n_agn_t1:
    RD = n.random.random(n_agn_t1)
    sel = (RD < N_qso * 1. / n_agn_t1)
    ids_qso = ids_type1[sel]


print("writing", time.time() - t0)
if status == 'create':
    grr = f.create_group('4MOST_COSMO_tracers')
    grr.create_dataset('ids_ELG', data=ids_elg)
    grr.create_dataset('ids_QSO', data=ids_qso)
    #grr.create_dataset('ids_LyAQSO', data = ids_qso )
    grr.create_dataset('ids_BG_lz', data=ids_lrg1)
    grr.create_dataset('ids_BG_hz', data=ids_lrg2)


f.close()
