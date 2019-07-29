"""
What it does
------------

Creates cluster images for each cluster in the light cone (healpix pixel of 13.4 deg2)

- 1. Create the number density + temperature profile: Eq. 7 + Table 3 of Ghirardini et al. 2019
- 2. interpolate scatter as a function of r/r500c using Table 2 for both. Correlate to coolness (I am not sure how to do that). What I have is a value between 0 and 1 that indicated the state of coolness. So we could multiply the scatter by a function of the coolness to obtain a smooth increase of scatter with 'disturbness'.
- 3. Normalize the \int_{0.01r500c}^{2r500c}{ne^2(r)} profile using the total luminosity that I with with the scaling relations (Bulbul et al. 2019) to create a X-ray surface brightness profile: SB(r)=Norm x ne^2(r)
- 4. Project on the sky: SB(theta), assuming a sphere. Convert directly T(r) to T(theta) with angular scale as a function of redshift.
- 5. Create an elliptical image (transformation that conserves the flux). I will take a set of discrete ellipticity values.
- 6. Slice it in elliptical shells and assign a temperature + Xspec APEC spectrum
- 7. create a base of decomposition of clusters into elliptical shells with X-ray spectra as a function of redshift to have less than 1000 images. By hand 10 redshifts x 5 shells x 4 ellipticities x 5 temperatures. With PCA (with non-negative decomposition).

References
----------

Ghirardini et al. 2019

Command to run
--------------

python3 004_6_cluster_images_simput.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

"""
import sys
import os
import time
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as p
import matplotlib
from scipy.integrate import quad
from scipy.interpolate import interp1d
from astropy_healpix import healpy
import astropy.io.fits as fits
import h5py
import numpy as n
print('CREATES SIMPUT CLUSTER FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()
#import astropy.io.fits as fits
#import healpy
# import all pathes
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})

env = sys.argv[1]  # 'MD04'
print(env)

write_spectra = False

# simulation setup
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

root_dir = os.path.join(os.environ[env])

dir_2_eRO_all = os.path.join(root_dir, "cat_eRO_CLU")
dir_2_SMPT = os.path.join(root_dir, "cat_CLU_SIMPUT")
dir_2_SMPT_image = os.path.join(root_dir, "cat_CLU_SIMPUT", 'cluster_images')
dir_2_SMPT_spectra = os.path.join(
    root_dir, "cat_CLU_SIMPUT", 'cluster_Xspectra')

if os.path.isdir(dir_2_SMPT_image) == False:
    os.system('mkdir -p ' + dir_2_SMPT_image)
if os.path.isdir(dir_2_SMPT_spectra) == False:
    os.system('mkdir -p ' + dir_2_SMPT_spectra)
if os.path.isdir(dir_2_SMPT) == False:
    os.system('mkdir -p ' + dir_2_SMPT)

fig_dir = os.path.join(
    os.environ['GIT_AGN_MOCK'],
    'figures',
    env,
    'clusters',
)
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)

# Ghirardini et al. 2019
# equation 7 f^2(x) = ne^2(x)
# x = R/R500c


def ne2(x, n0, rc, alpha, beta, rs, epsilon, gamma=3.): return n0**2. * ((x / rc)**(- alpha) / (1 +
                                                                                                x**2. / rc**2.)**(3. * beta - alpha / 2.)) * (1. / (1. + x**gamma / rs**gamma) ** (epsilon / gamma))


# cool core profile
n0_cc = 10**-3.9
rc_cc = 10**-3.2
rs_cc = 10**0.17
alpha_cc = 0.80
beta_cc = 0.49
epsilon_cc = 4.67


def ne2_cc(x): return ne2(
    x,
    n0_cc,
    rc_cc,
    alpha_cc,
    beta_cc,
    rs_cc,
    epsilon_cc)


# Non cool core profile
n0_ncc = 10**-4.9
rc_ncc = 10**-2.7
rs_ncc = 10**-0.51
alpha_ncc = 0.70
beta_ncc = 0.39
epsilon_ncc = 2.6


def ne2_ncc(x): return ne2(
    x,
    n0_ncc,
    rc_ncc,
    alpha_ncc,
    beta_ncc,
    rs_ncc,
    epsilon_ncc)


# surface brightness integration
def fun_cc(x, xi=0.1): return x * ne2_cc(x) * (x**2. - xi**2.)**(-0.5)


def fun_ncc(x, xi=0.1): return x * ne2_ncc(x) * (x**2. - xi**2.)**(-0.5)


def itg_cc(xi): return xi * quad(fun_cc, xi, 2, args=(xi))[0]


def itg_ncc(xi): return xi * quad(fun_ncc, xi, 2, args=(xi))[0]


# angular bins, in fraction of theta_500c
# n.array([0.001, 0.01, 0.07, 0.13, 0.21, 0.31, 0.46, 0.72, 1.15])
xis = 10**n.arange(-3, 0.3, 0.1)
#delta_xis = xis[1:]-xis[:-1]
F_cc_xi = n.array([itg_cc(xi_i) for xi_i in xis])
F_ncc_xi = n.array([itg_ncc(xi_i) for xi_i in xis])
# interpolation, x axis: 1 rvir = 40 pixels
#                               = 40 * 5.5e-04 # degrees
rc_pixels = 40.
rc_arcseconds = rc_pixels * 5.5e-04 * 3600  # arc seconds
F_cc_xi_itp = interp1d(
    n.hstack((0., xis, 10.)) * 40,
    n.hstack((F_cc_xi[0], F_cc_xi, 0.))
)
F_ncc_xi_itp = interp1d(
    n.hstack((0., xis, 10.)) * 40,
    n.hstack((F_ncc_xi[0], F_ncc_xi, 0.))
)


#fig_out = os.path.join(fig_dir, 'SB_profiles.png')
#p.figure(1, (6.,6.))
# p.axes([0.18,0.15,0.75,0.78])
#p.plot(xis, F_cc_xi, label='cc')
#p.plot(xis, F_ncc_xi, label='ncc')
# p.legend(frameon=False)
# p.xlabel(r'$\theta/\theta_{500c}$')
#p.ylabel('surface brightness')
# p.grid()
# p.ylim((0,1.1))
# p.yscale('log')
# p.xscale('log')
# p.savefig(fig_out)
# p.clf()


# Create a circular image
# normalize the flux to the catalogs
# truncate to 2xr500c
#

n_pixel = 120
matrix = n.zeros((n_pixel, n_pixel))
xxx = (n.arange(n_pixel) - n_pixel / 2.)  # arc seconds
x_matrix, y_matrix = n.meshgrid(xxx, xxx)


def write_img(matrix, name='spherical_cc'):
    prihdr = fits.Header()
    prihdr['HDUCLASS'] = 'HEASARC/SIMPUT'
    prihdr['HDUCLAS1'] = 'IMAGE'
    prihdr['HDUVERS'] = '1.1.0'
    prihdr['EXTNAME'] = 'IMAGE'
    prihdr['CTYPE1'] = ('RA---TAN', 'first axis (column) is Right Ascension')
    prihdr['CRPIX1'] = (n_pixel / 2., 'middle pixel of array in col direction')
    prihdr['CRVAL1'] = (0, 'Dec of this middle pixel, in degrees')
    prihdr['CDELT1'] = (-5.5e-04,
                        'move 1column forward,decrease RA by CDELT1/deg')
    prihdr['CROTA1'] = 0
    prihdr['CUNIT1'] = 'deg'
    prihdr['CTYPE2'] = ('DEC--TAN', 'first axis (column) is Declination')
    prihdr['CRPIX2'] = (n_pixel / 2., 'middle pixel of array in row direction')
    prihdr['CRVAL2'] = (0, 'RA of this middle pixel, in degrees')
    prihdr['CDELT2'] = (
        5.5e-04,
        'move 1column forward,increase Dec by CDELT1/deg')
    prihdr['CROTA2'] = 0
    prihdr['CUNIT2'] = 'deg'
    prihdr['EQUINOX'] = 2000
    prihdr['RADECSYS'] = 'FK5'
    prihdu = fits.PrimaryHDU(matrix, header=prihdr)
    #out = os.path.join(os.environ['MD10'], 'sixte', 'images', name+'.fits')
    out = os.path.join(dir_2_SMPT_image, name + '.fits')
    if os.path.isfile(out):
        os.remove(out)
    prihdu.writeto(out, clobber=True)


if write_spectra:
    # make them elliptical
    # Create an elliptical image, 4 ellipticities x 2 profiles ?

    # b_a=1.0
    r_matrix = ((x_matrix)**2 + y_matrix**2)**0.5
    cc_matrix = F_cc_xi_itp(r_matrix)
    ncc_matrix = F_ncc_xi_itp(r_matrix)
    cc_matrix[r_matrix > 60] = 0
    ncc_matrix[r_matrix > 60] = 0
    cc_matrix = cc_matrix / n.sum(cc_matrix)
    ncc_matrix = ncc_matrix / n.sum(ncc_matrix)
    write_img(cc_matrix, name='elliptical_ba_1p00_cc')
    write_img(ncc_matrix, name='elliptical_ba_1p00_ncc')

    # b_a=0.71
    r_matrix = ((x_matrix / 0.7)**2 + y_matrix**2)**0.5
    cc_matrix = F_cc_xi_itp(r_matrix)
    ncc_matrix = F_ncc_xi_itp(r_matrix)
    cc_matrix[r_matrix > 60] = 0
    ncc_matrix[r_matrix > 60] = 0
    cc_matrix = cc_matrix / n.sum(cc_matrix)
    ncc_matrix = ncc_matrix / n.sum(ncc_matrix)
    write_img(cc_matrix, name='elliptical_ba_0p70_cc')
    write_img(ncc_matrix, name='elliptical_ba_0p70_ncc')

    # b_a=0.5
    r_matrix = ((x_matrix / 0.5)**2 + y_matrix**2)**0.5
    cc_matrix = F_cc_xi_itp(r_matrix)
    ncc_matrix = F_ncc_xi_itp(r_matrix)
    cc_matrix[r_matrix > 60] = 0
    ncc_matrix[r_matrix > 60] = 0
    cc_matrix = cc_matrix / n.sum(cc_matrix)
    ncc_matrix = ncc_matrix / n.sum(ncc_matrix)
    write_img(cc_matrix, name='elliptical_ba_0p50_cc')
    write_img(ncc_matrix, name='elliptical_ba_0p50_ncc')

    # b_a=0.25
    r_matrix = ((x_matrix / 0.25)**2 + y_matrix**2)**0.5
    cc_matrix = F_cc_xi_itp(r_matrix)
    ncc_matrix = F_ncc_xi_itp(r_matrix)
    cc_matrix[r_matrix > 60] = 0
    ncc_matrix[r_matrix > 60] = 0
    cc_matrix = cc_matrix / n.sum(cc_matrix)
    ncc_matrix = ncc_matrix / n.sum(ncc_matrix)
    write_img(cc_matrix, name='elliptical_ba_0p25_cc')
    write_img(ncc_matrix, name='elliptical_ba_0p25_ncc')

# link to names in simput file.
N_pixels = healpy.nside2npix(8)
for HEALPIX_id in n.arange(N_pixels):
    #HEALPIX_id = 358
    path_2_eRO_all_catalog = os.path.join(
        dir_2_eRO_all, str(HEALPIX_id).zfill(6) + '.fit')
    path_2_SMPT_catalog = os.path.join(
        dir_2_SMPT, 'SIMPUT_' + str(HEALPIX_id).zfill(6) + '.fit')

    hd_all = fits.open(path_2_eRO_all_catalog)
    N_clu_all = len(hd_all[1].data['ra'])
    ra_array = hd_all[1].data['ra']
    dec_array = hd_all[1].data['dec']
    redshift = hd_all[1].data['redshift_R']
    coolness = hd_all[1].data['coolness']
    HALO_M500c = hd_all[1].data['HALO_M500c']
    HALO_b_to_a_500c = hd_all[1].data['HALO_b_to_a_500c']
    HALO_c_to_a_500c = hd_all[1].data['HALO_c_to_a_500c']
    FX_soft_attenuated = hd_all[1].data['FX_soft_attenuated']
    HALO_rs = hd_all[1].data['HALO_rs']
    HALO_rvir = hd_all[1].data['HALO_rvir']
    HALO_Mvir = hd_all[1].data['HALO_Mvir']
    kT = hd_all[1].data['TX_cin']

    rd_all = n.random.rand(N_clu_all)
    orientation = n.random.rand(N_clu_all) * 180.  # IMGROTA

    arcsec_p_kpc = cosmo.arcsec_per_kpc_comoving(redshift).value
    pixel_rescaling = rc_arcseconds / (arcsec_p_kpc * HALO_rvir)  # IMGSCALE
    #print('pixel_rescaling', pixel_rescaling, pixel_rescaling.min(), pixel_rescaling.max())
    pixel_rescaling[pixel_rescaling < 0.001] = 0.001
    # NOW ASSIGNS TEMPLATES BASED ON THE HALO PROPERTIES
    ba_025 = (HALO_b_to_a_500c < 0.25)
    ba_050 = (HALO_b_to_a_500c >= 0.25) & (HALO_b_to_a_500c < 0.5)
    ba_075 = (HALO_b_to_a_500c >= 0.50) & (HALO_b_to_a_500c < 0.75)
    ba_100 = (HALO_b_to_a_500c >= 0.75)

    # https://arxiv.org/abs/1703.08690 Felipe Andrade-Santos
    # fraction of cool core is 44% pm 7% up to redshift 0.35
    coolcore = (coolness > 0.56)  # relaxed
    non_coolcore = (coolness <= 0.56)  # non-relaxed

    template = n.zeros(N_clu_all).astype('U100')
    template[template ==
             "0.0"] = "cluster_images/elliptical_ba_0p25_cc.fits[SPECTRUM][#row==1]"

    template[ba_025 & coolcore] = """cluster_images/elliptical_ba_0p25_cc.fits[IMAGE]"""
    template[ba_050 & coolcore] = """cluster_images/elliptical_ba_0p50_cc.fits[IMAGE]"""
    template[ba_075 & coolcore] = """cluster_images/elliptical_ba_0p70_cc.fits[IMAGE]"""
    template[ba_100 & coolcore] = """cluster_images/elliptical_ba_1p00_cc.fits[IMAGE]"""
    template[ba_025 & non_coolcore] = """cluster_images/elliptical_ba_0p25_ncc.fits[IMAGE]"""
    template[ba_050 & non_coolcore] = """cluster_images/elliptical_ba_0p50_ncc.fits[IMAGE]"""
    template[ba_075 & non_coolcore] = """cluster_images/elliptical_ba_0p70_ncc.fits[IMAGE]"""
    template[ba_100 & non_coolcore] = """cluster_images/elliptical_ba_1p00_ncc.fits[IMAGE]"""

    # NOW links to the grid of SPECTRA

    kt_arr = n.array([0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 10.0])
    z_arr = n.hstack((n.array([0., 0.05]), n.arange(0.1, 1.6, 0.1)))

    indexes_kt = n.array([(n.abs(kT_val - kt_arr)).argmin()
                          for kT_val in 10**kT])
    kT_values = kt_arr[indexes_kt]

    indexes_z = n.array([(n.abs(z_val - z_arr)).argmin()
                         for z_val in redshift])
    z_values = z_arr[indexes_z]

    def tpl_name(temperature, redshift): return 'cluster_Xspectra/cluster_spectrum_10kT_' + str(int(
        temperature * 10)).zfill(4) + '_100z_' + str(int(redshift * 100)).zfill(4) + '.fits[SPECTRUM][#row==1]'

    spec_names = n.zeros(N_clu_all).astype('U200')
    # spec_names[template=="0.0"] =
    # "cluster_Xspectra/cluster_spectrum_10kT_0100_100z_0150.fits[SPECTRUM][#row==1]"

    for jj, (kT_values_ii, z_values_ii) in enumerate(zip(kT_values, z_values)):
        spec_names[jj] = tpl_name(kT_values_ii, z_values_ii)

    hdu_cols = fits.ColDefs([
        # ,fits.Column(name="SRC_NAME" , format='10A', unit='', array =  n.arange(len(hd_all[1].data['ra'])).astype('str') )
        fits.Column(name="SRC_ID", format='K', unit='', array=(n.arange(N_clu_all) + 4e8).astype('int')), fits.Column(name="RA", format='D', unit='deg', array=ra_array), fits.Column(name="DEC", format='D', unit='deg', array=dec_array), fits.Column(name="E_MIN", format='D', unit='keV', array=n.ones(N_clu_all) * 0.5), fits.Column(name="E_MAX", format='D', unit='keV', array=n.ones(N_clu_all) * 2.0), fits.Column(name="FLUX", format='D', unit='erg/s/cm**2', array=FX_soft_attenuated), fits.Column(name="IMAGE", format='100A', unit='', array=template), fits.Column(name="SPECTRUM", format='100A', unit='', array=spec_names)        # ,fits.Column(name= "n_energy_bins"  , format='K', unit='', array = data_n_e_b)
        , fits.Column(name="IMGROTA", format='K', unit='deg', array=orientation), fits.Column(name="IMGSCAL", format='K', unit='', array=pixel_rescaling)
    ])

    hdu = fits.BinTableHDU.from_columns(hdu_cols)

    hdu.name = 'SRC_CAT'
    hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
    hdu.header['HDUCLAS1'] = 'SRC_CAT'
    hdu.header['HDUVERS'] = '1.1.0'
    hdu.header['RADESYS'] = 'FK5'
    hdu.header['EQUINOX'] = 2000.0

    outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])
    if os.path.isfile(path_2_SMPT_catalog):
        os.system("rm " + path_2_SMPT_catalog)
    outf.writeto(path_2_SMPT_catalog, overwrite=True)
    print(path_2_SMPT_catalog, 'written', time.time() - t0)
