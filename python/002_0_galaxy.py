"""
What it does
------------

Create a file with galaxy properties.

Computes
 - stellar mass following Moster et al. 2013, 2018
 - star formation rates following the main sequence of galaxies following Whitaker et al 2012
 - star formation rate for quiescent galaxies following fits on COSMOS data, see Comparat et al. in prep
 - X-ray luminosities of galaxies after Aird et al. 2018
 - r-band magnitude in the SDSS AB system doing abundance matching with the luminosity function of Loveday et al. 2016

References
----------

Whitaker 2012, 2014
https://ui.adsabs.harvard.edu/abs/2012ApJ...754L..29W
https://ui.adsabs.harvard.edu/abs/2014ApJ...795..104W
Loveday et al. 2016

Ilber et al. 2013


Command to run
--------------

python3 002_0_galaxy.py environmentVAR fileBasename

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. MD10
It will then work in the directory : $environmentVAR/hlists/fits/

fileBasename: base name of the file e.g. all_1.00000

It will find the input files :
$environmentVAR/hlists/fits/${fileBasename}.fits
$environmentVAR/hlists/fits/${fileBasename}_coordinates.h5

And will write outputs in h5 format:
$environmentVAR/hlists/fits/${fileBasename}_galaxy.h5

Figures (Optional)
if the variable 'make_figure' is set to True, then figures will be created in the git repo here :
$GIT_AGN_MOCK/figures/$environmentVAR/galaxy/
It makes relevant histograms and scatter plots for all columns present in the new file.

"""
import sys
import os
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp1d
import h5py
import astropy.io.fits as fits
import numpy as n
import time
print('Adds galaxy properties ')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


env = sys.argv[1]  # 'MD04'
baseName = sys.argv[2]  # "sat_0.62840"
print(env, baseName)
make_figure = True
make_figure = False

# import all pathes
test_dir = os.path.join(os.environ[env], 'hlists', 'fits')

path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
path_2_coordinate_file = os.path.join(test_dir, baseName + '_coordinates.h5')
path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.h5')

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

print('opens light cone ')
f1 = fits.open(path_2_light_cone)
Mvir = f1[1].data['Mvir'] / h
log_vmax = n.log10(f1[1].data['vmax'])
N_obj = len(Mvir)
f1.close()

print('opens coordinates ')
f2 = h5py.File(path_2_coordinate_file, 'r')
zz = f2['/coordinates/redshift_R'][:]
N_halos = len(zz)
f2.close()

N_galaxies = len(zz)
volume = (
    cosmo.comoving_volume(
        n.max(zz)) -
    cosmo.comoving_volume(
        n.min(zz))).value
print(volume)

# STELLAR MASS
# Equations 1 of Comparat et al. 2019


def meanSM(Mh, z): return n.log10(Mh * 2. * (0.0351 - 0.0247 * z / (1. + z)) / ((Mh / (10**(11.79 + 1.5 * z / (1. + z))))
                                                                                ** (- 0.9 + 0.5 * z / (1. + z)) + (Mh / (10**(11.79 + 1.5 * z / (1. + z))))**(0.67 + 0.2 * z / (1. + z))))


mean_SM = meanSM(Mvir, zz)
rds = norm.rvs(loc=0, scale=0.15, size=N_halos)
mass = mean_SM + rds  # fun(mean_SM)
print('masses', mass[:10], time.time() - t0)

# STAR FORMATION RATE (valid only for star forming galaxies !)
sfr = n.zeros_like(zz)
# whitaker et al 2012, Eq. 1,2,3.
mean_SFR = (0.70 - 0.13 * zz) * (mass - 10.5) + \
    0.38 + 1.14 * zz - 0.19 * zz**2.
rds2 = norm.rvs(loc=0, scale=0.34, size=N_halos)
log_sfr = mean_SFR + rds2
print('log sfr', log_sfr[:10], time.time() - t0)

# quiescent fraction fitted on COSMOS, Ilbert et al. 2013 v2 catalogue


def scatter_Qf(z): return - 0.45 * (1 + z) + 1.54


def log10_M0_Qf(z): return 9.71 + 0.78 * (1 + z)


def fraction_Qf(mass, z): return 0.5 + 0.5 * \
    erf((mass - log10_M0_Qf(z)) / scatter_Qf(z))


rds_Qf = n.random.random(N_obj)

frac = fraction_Qf(mass, zz)
SF = (rds_Qf > frac)
QU = (SF == False)

# mass-SFR sequence for the quiescent
sfr_Q = n.zeros_like(zz)


def beta_z(z): return -0.57 * z + 1.43


def alpha_z(z): return 6.32 * z - 16.26


def mean_SFR_Q(mass, z): return mass * beta_z(z) + alpha_z(z)


def scale_z(z): return -0.34 * z + 0.99


rds2 = norm.rvs(loc=0, scale=1., size=len(zz[QU])) * scale_z(zz[QU])

log_sfr_Q = mean_SFR_Q(mass[QU], zz[QU]) + rds2
# change SFR for the quiesent selection
log_sfr[QU] = log_sfr_Q

print('N QU, log SFR[:10]', len(log_sfr[QU]),
      log_sfr[QU][:10], time.time() - t0)
print('N SF, log SFR[:10]', len(log_sfr[SF]),
      log_sfr[SF][:10], time.time() - t0)

# Hard X-ray emission, after Aird et al. 2018


def galaxy_lx(redshift, mass, sfr):
    return 10**(28.81) * (1 + redshift)**(3.9) * mass + \
        10**(39.5) * (1 + redshift)**(0.67) * sfr**(0.86)


gal_LX = galaxy_lx(zz, 10**mass, 10**log_sfr)
print('gal_LX', gal_LX[:10])

# galaxy LF abundance matching to Loveday et al. 2016
ngal, M, phi, Err = n.loadtxt(os.path.join(
    os.environ['GIT_AGN_MOCK'], 'data', 'LF_loveday_2015', 'lf.txt'), unpack=True)
Schechter_M_z = interp1d(M, phi * 0.7**3.)
mrange = n.arange(-24.6, -12.2, 0.01)
total_number = n.cumsum(Schechter_M_z(mrange)) * volume
print('total_number', total_number)

rds = norm.rvs(loc=0, scale=0.15, size=N_galaxies)
print('N_galaxies', N_galaxies)
vmax_galaxy = 10**(log_vmax + rds)

x_itp = n.hstack((0.,
                  total_number[total_number > 0.1],
                  total_number[total_number > 0.1][-1] * 10))
y_itp = n.hstack((mrange[total_number > 0.1][0],
                  mrange[total_number > 0.1],
                  mrange[total_number > 0.1][-1]))
itp_mags = interp1d(x_itp, y_itp)
mags_2_assign = itp_mags(n.arange(N_galaxies) + 1)
print('mags_2_assign', mags_2_assign, mags_2_assign.shape)
id_sort_vmax = n.argsort(vmax_galaxy)[::-1]  # id=0 the lowest vmax
mag_abs_r = mags_2_assign[id_sort_vmax]
print('mag_abs_r mag_r', mag_abs_r[:10], time.time() - t0)

# compute distance modulus
z_array = n.arange(n.min(zz) - 0.05, n.max(zz) + 0.05, 0.0001)
dist_mods = interp1d(z_array, cosmo.distmod(z_array).value)
mag_r = mag_abs_r + dist_mods(zz)
print('mag_r', mag_r[:10], time.time() - t0)

"""
# mass size relation

# https://arxiv.org/pdf/1411.6355.pdf
# Table 2 and 3. r mag line for sersic selection
# equation 2

# http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
b4 = 7.669
b1 = 1.678
f_14_dev = lambda r12 : gammainc(8, b4*(0.7/r12)**(1./6.))
f_14_exp = lambda r12 : gammainc(2, b1*(0.7/r12)**(1./2.)) # From Raichoor 2017
f_14_test = lambda r12, nn : gammainc(2, b1*(0.7/r12)**(1./nn))

re_dev = lambda M_star : 0.16*(M_star)**(0.1) * (1+M_star/(2.42*10**(10)))**(0.76-0.1)
re_exp = lambda M_star : 0.08*(M_star)**(0.16) * (1+M_star/(17.1*10**(10)))**(0.81-0.16)
# interpolates projected size on the sky
radius_2_arcsec = interp1d(n.arange(0.00001,6.5,0.001), lc_setup.cosmoMD.arcsec_per_kpc_proper( n.arange(0.00001,6.5,0.001) ).value)

radius_dev = re_dev(stellar_mass)* radius_2_arcsec(f1['/sky_position/redshift_R'].value[gal][ok])
radius_exp = re_exp(stellar_mass)* radius_2_arcsec(f1['/sky_position/redshift_R'].value[gal][ok])

radius = radius_exp

dm_dev = -2.5*n.log10(f_14_dev(radius))
dm_exp = -2.5*n.log10(f_14_exp(radius))
# initialize the fibermag with the exponential profile
fiber_mag = rmag + dm_exp
"""

f = h5py.File(path_2_galaxy_file, "a")
f.attrs['file_name'] = os.path.basename(path_2_galaxy_file)
f.attrs['creator'] = 'JC'

# writes the results
halo_data = f.create_group('galaxy')

ds = halo_data.create_dataset('SMHMR_mass', data=mass)
ds.attrs['units'] = 'log10(stellar mass/[Msun])'

ds = halo_data.create_dataset('star_formation_rate', data=log_sfr)
ds.attrs['units'] = 'log10(SFR/[Msun/yr])'

ds = halo_data.create_dataset('is_quiescent', data=QU)
ds.attrs['units'] = 'boolean'

ds = halo_data.create_dataset('LX_hard', data=n.log10(gal_LX))
ds.attrs['units'] = 'log10(L_X/[2-10keV, erg/s])'

ds = halo_data.create_dataset('mag_abs_r', data=mag_abs_r)
ds = halo_data.create_dataset('mag_r', data=mag_r)
ds.attrs['units'] = 'mag AB'

f.close()

### Option: FIGURES ###
if make_figure:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': 14})
    import matplotlib.pyplot as p

    fig_dir = os.path.join(
        os.environ['GIT_AGN_MOCK'],
        'figures',
        env,
        'galaxy',
    )
    if os.path.isdir(fig_dir) == False:
        os.system('mkdir -p ' + fig_dir)

    h5_file = os.path.join(path_2_coordinate_file)
    f = h5py.File(h5_file, "r")
    zr = f['/coordinates/redshift_R'][:]
    f.close()

    N_obj = len(mass)
    sel = (rds < 1e5 / N_obj)

    # histograms

    fig_out = os.path.join(fig_dir, 'SMHMR_mass_hist_' + baseName + '.png')

    p.figure(1, (6., 5.5))
    p.tight_layout()
    X = mass
    p.hist(X, histtype='step', label='all', rasterized=True, lw=4)
    X = mass[QU]
    p.hist(X, histtype='step', label='QU', rasterized=True, lw=2, ls='dashed')
    X = mass[SF]
    p.hist(X, histtype='step', label='SF', rasterized=True, lw=2)
    p.title(baseName)
    p.xlabel('SMHMR_mass')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.legend(frameon=False, loc=0)
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(
        fig_dir,
        'star_formation_rate_hist_' +
        baseName +
        '.png')

    p.figure(1, (6., 5.5))
    p.tight_layout()
    X = log_sfr
    p.hist(X, histtype='step', label='all', rasterized=True, lw=4)
    X = log_sfr[SF]
    p.hist(X, histtype='step', label='SF', rasterized=True, lw=2, ls='dashed')
    X = log_sfr[QU]
    p.hist(X, histtype='step', label='QU', rasterized=True, lw=2)
    p.title(baseName)
    p.xlabel('galaxy_star_formation_rate')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.legend(frameon=False, loc=0)
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'LX_hard_hist_' + baseName + '.png')

    X = n.log10(gal_LX)

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('galaxy_LX_hard')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'mag_abs_r_hist_' + baseName + '.png')

    X = mag_abs_r

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('mag_abs_r')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'mag_r_hist_' + baseName + '.png')

    X = mag_r

    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.hist(X, histtype='step', rasterized=True, lw=4)
    p.title(baseName)
    p.xlabel('galaxy_mag_r')
    p.ylabel('Counts')
    p.grid()
    p.yscale('log')
    p.savefig(fig_out)
    p.clf()

    # 2D plots

    fig_out = os.path.join(fig_dir, 'mag_r_vs_zr_' + baseName + '.png')

    X = zr[sel]
    Y = mag_r[sel]
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k,', rasterized=True)
    p.title(baseName)
    p.xlabel('redshift R')
    p.ylabel('mag_r')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'mag_abs_r_vs_zr_' + baseName + '.png')

    X = zr[sel]
    Y = mag_abs_r[sel]
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k,', rasterized=True)
    p.title(baseName)
    p.xlabel('redshift R')
    p.ylabel('mag_abs_r')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(fig_dir, 'SMHMR_mass_vs_zr_' + baseName + '.png')

    X = zr[sel]
    Y = mass[sel]
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k,', rasterized=True)
    p.title(baseName)
    p.xlabel('redshift R')
    p.ylabel('log10 SMHMR_mass')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(
        fig_dir,
        'star_formation_rate_vs_zr_' +
        baseName +
        '.png')

    X = zr[sel]
    Y = log_sfr[sel]
    p.figure(1, (6., 5.5))
    p.tight_layout()
    p.plot(X, Y, 'k+', rasterized=True)
    p.title(baseName)
    p.xlabel('redshift R')
    p.ylabel('log10 star_formation_rate')
    p.grid()
    p.savefig(fig_out)
    p.clf()

    fig_out = os.path.join(
        fig_dir,
        'star_formation_rate_vs_mass_' +
        baseName +
        '.png')

    p.figure(1, (6., 5.5))
    p.tight_layout()
    X = mass[sel]
    Y = log_sfr[sel]
    p.plot(X, Y, 'b+', rasterized=True)
    X = mass[sel & QU]
    Y = log_sfr[sel & QU]
    p.plot(X, Y, 'r+', label='QU', rasterized=True)
    p.title(baseName)
    p.legend(frameon=False, loc=0)
    p.xlabel('log10(mass)')
    p.ylabel('log10 star_formation_rate')
    p.grid()
    p.savefig(fig_out)
    p.clf()
