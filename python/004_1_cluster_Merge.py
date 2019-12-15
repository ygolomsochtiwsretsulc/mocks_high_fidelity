"""
What it does
------------

Merges the cluster catalogs into a single full sky catalog in fits format

References
----------

Command to run
--------------

python3 004_1_cluster_Merge.py environmentVAR laptop

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/


laptop: if running on the server or on the laptop. Redirects properly to the stilts command

Dependencies
------------

topcat/stilts
import time, os, sys, numpy, scipy, astropy, h5py, astropy_healpix, matplotlib



"""
import glob
from astropy_healpix import healpy
import sys
import os
import time
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.constants as cc
import astropy.io.fits as fits
from scipy.special import erf
from scipy.stats import norm
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import h5py
import numpy as n
print('Creates Cluster file ')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()


# import all pathes

env = sys.argv[1]  # 'MD04'
laptop = sys.argv[2]  # 'True'

print(env)
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

if laptop == "True":
    stilts_cmd = 'java -jar /home/comparat/software/stilts.jar'
else:
    stilts_cmd = 'stilts'

# initializes pathes to files
test_dir = os.path.join(os.environ[env], 'fits')
tmp_dir = os.path.join(os.environ[env], 'fits', 'tmp')
if os.path.isdir(tmp_dir) == False:
    os.system('mkdir -p ' + tmp_dir)

path_2_CLU_catalog = os.path.join(test_dir, env + '_eRO_CLU.fit')

if os.path.isfile(path_2_CLU_catalog):
    os.system("rm " + path_2_CLU_catalog)


path_2_CLU_files = sorted(
    n.array(
        glob.glob(
            os.path.join(
                test_dir,
                '*_CLU.fits'))))

baseNames = n.array([os.path.basename(path_2_CLU_file)[:-7]
                     for path_2_CLU_file in path_2_CLU_files])


def get_data(
        path_2_light_cone,
        path_2_coordinate_file,
        path_2_galaxy_file,
        path_2_CLU_file):
    f3 = h5py.File(path_2_CLU_file, 'r')
    clu = f3['/CLUSTERS/ids_cluster'].value
    LX_soft_cin = f3['/CLUSTERS/LX_soft_cin'].value
    LX_soft_cex = f3['/CLUSTERS/LX_soft_cex'].value
    TX_cin = f3['/CLUSTERS/TX_cin'].value
    TX_cex = f3['/CLUSTERS/TX_cex'].value
    MICM_cex = f3['/CLUSTERS/MICM_cex'].value
    FX_soft = f3['/CLUSTERS/FX_soft'].value
    FX_soft_attenuated = f3['/CLUSTERS/FX_soft_attenuated'].value
    detected = f3['/CLUSTERS/detected'].value
    scatter_1 = f3['/CLUSTERS/scatter_1'].value
    coolness = f3['/CLUSTERS/coolness'].value
    f3.close()

    f1 = h5py.File(path_2_galaxy_file, 'r')
    galaxy_stellar_mass = f1['galaxy/SMHMR_mass'].value[clu]
    galaxy_star_formation_rate = f1['galaxy/star_formation_rate'].value[clu]
    galaxy_LX_hard = f1['galaxy/LX_hard'].value[clu]
    galaxy_mag_r = f1['galaxy/mag_r'].value[clu]
    galaxy_mag_abs_r = f1['galaxy/mag_abs_r'].value[clu]
    f1.close()

    f2 = h5py.File(path_2_coordinate_file, 'r')
    ra = f2[1].data['ra'].value[clu]
    dec = f2[1].data['dec'].value[clu]
    zzr = f2[1].data['redshift_R'].value[clu]
    zzs = f2[1].data['redshift_S'].value[clu]
    dL_cm = f2[1].data['dL'].value[clu]
    galactic_NH = f2['nH'].value[clu]
    galactic_ebv = f2[1].data['ebv'].value[clu]
    g_lat = f2[1].data['g_lat'].value[clu]
    g_lon = f2[1].data['g_lon'].value[clu]
    ecl_lat = f2[1].data['ecl_lat'].value[clu]
    ecl_lon = f2[1].data['ecl_lon'].value[clu]
    N_galaxies = len(zzr)
    f2.close()

    HALO_cat = fits.open(path_2_light_cone)[1].data[clu]

    hdu_cols = fits.ColDefs([
        ##
        # Coordinates
        ##
        fits.Column(
            name="ra", format='D', unit='degree', array=ra), fits.Column(
            name="dec", format='D', unit='degree', array=dec), fits.Column(
            name="g_lat", format='D', unit='degree', array=g_lat), fits.Column(
                name="g_lon", format='D', unit='degree', array=g_lon), fits.Column(
                    name="ecl_lat", format='D', unit='degree', array=ecl_lat), fits.Column(
                        name="ecl_lon", format='D', unit='degree', array=ecl_lon)        # distances
        , fits.Column(name="redshift_R", format='D', unit='real space', array=zzr), fits.Column(name="redshift_S", format='D', unit='redshift space', array=zzs), fits.Column(name="dL_cm", format='D', unit='cm', array=dL_cm)        # extinction maps
        # Galaxy properties
        # AGN properties
        # ,fits.Column(name= "AGN_LX_soft"          , format='D', unit='Luminosity/[erg/s] 0.5-2 keV', array = agn_in_clu_LX_soft  )
        , fits.Column(name="galactic_NH", format='D', unit='cm-2', array=galactic_NH), fits.Column(name="galactic_ebv", format='D', unit='mag', array=galactic_ebv), fits.Column(name="galaxy_stellar_mass", format='D', unit='log10(M/[M_sun])', array=galaxy_stellar_mass), fits.Column(name="galaxy_star_formation_rate", format='D', unit='log10(SFR/[M_sun/year])', array=galaxy_star_formation_rate), fits.Column(name="galaxy_LX_hard", format='D', unit='log10(LX (2-10keV)/[erg/s])', array=galaxy_LX_hard), fits.Column(name="galaxy_mag_r", format='D', unit='mag', array=galaxy_mag_r), fits.Column(name="galaxy_mag_abs_r", format='D', unit='mag err', array=galaxy_mag_abs_r)        # ,fits.Column(name= "AGN_FX_soft"          , format='D', unit='Flux/[erg/cm2/s] 0.5-2 keV',   array = agn_in_clu_FX_soft_attenuated    )
        # ,fits.Column(name= "AGN_LX_hard"          , format='D', unit='Luminosity/[erg/s] 2-10 keV',  array = agn_in_clu_LX_hard   )
        # ,fits.Column(name= "AGN_FX_hard"          , format='D', unit='Flux/[erg/cm2/s] 2-10 keV',    array = agn_in_clu_FX_hard     )
        # ,fits.Column(name= "AGN_SDSS_r_magnitude" , format='D', unit='mag',                          array = agn_in_clu_SDSS_r_AB_attenuated          )
        # ,fits.Column(name= "AGN_Nh"               , format='D', unit='log10(Nh/[cm-2])',             array = agn_in_clu_logNH               )
        # ,fits.Column(name= "AGN_random_number"    , format='D', unit='',                             array = agn_in_clu_random       )
        # ,fits.Column(name= "AGN_type"             , format='D', unit='X/opt type: 11, 12, 21, 22',   array = agn_in_clu_agn_type           )
        ##
        # Dark matter halo
        ##
        , fits.Column(name="HALO_M200c", format='D', unit='log10(M/[M_sun])', array=HALO_cat["M200c"] / h), fits.Column(name="HALO_M500c", format='D', unit='log10(M/[M_sun])', array=HALO_cat["M500c"] / h), fits.Column(name="HALO_Mvir", format='D', unit='log10(M/[M_sun])', array=HALO_cat["Mvir"] / h), fits.Column(name="HALO_Acc_Rate_1Tdyn", format='D', unit='Msun/yr', array=HALO_cat["Acc_Rate_1_Tdyn"]), fits.Column(name="HALO_b_to_a_500c", format='D', unit='axis ratio', array=HALO_cat["b_to_a_500c"]), fits.Column(name="HALO_c_to_a_500c", format='D', unit='axis ratio', array=HALO_cat["c_to_a_500c"]), fits.Column(name="HALO_rs", format='D', unit='kpc', array=HALO_cat["rs"]), fits.Column(name="HALO_rvir", format='D', unit='kpc', array=HALO_cat["rvir"]), fits.Column(name="HALO_vmax", format='D', unit='km/s', array=HALO_cat["vmax"]), fits.Column(name="x", format='D', unit='Mpc', array=HALO_cat["x"]), fits.Column(name="y", format='D', unit='Mpc', array=HALO_cat["y"]), fits.Column(name="z", format='D', unit='Mpc', array=HALO_cat["z"])        #
        , fits.Column(name="halo_id", format='K', unit='', array=HALO_cat['id'])        # ,fits.Column(name=  "pid"      , format='K', unit='', array = HALO_cat['pid'] )
        ##
        # cluster properties
        ##
        # ,fits.Column(name='CLU_FX_soft',            unit='erg/cm2/s', format='D', array= FX_soft  )
        # ,fits.Column(name='CLU_FX_soft_attenuated', unit='erg/cm2/s', format='D', array= CLU_FX_soft_attenuated  )
        # ,fits.Column(name='CLU_LX_soft',            unit='erg/s', format='D',  array=CLU_LX_soft  )
        , fits.Column(name='LX_soft_cin', unit='log10(L_X/[0.5-2keV, cin, erg/s])', format='D', array=LX_soft_cin), fits.Column(name='LX_soft_cex', unit='log10(L_X/[0.5-2keV, cex, erg/s])', format='D', array=LX_soft_cex), fits.Column(name='TX_cin', unit='log10(T cin [keV])', format='D', array=TX_cin), fits.Column(name='TX_cex', unit='log10(T cex [keV])', format='D', array=TX_cex), fits.Column(name='MICM_cex', unit='log10(M_ICM [Msun]))', format='D', array=MICM_cex), fits.Column(name='FX_soft', unit='F_X / [0.5-2keV, erg/cm2/s]', format='D', array=FX_soft), fits.Column(name='FX_soft_attenuated', unit='F_X / [0.5-2keV, erg/cm2/s]', format='D', array=FX_soft_attenuated), fits.Column(name='detected', unit='boolean', format='L', array=detected), fits.Column(name='scatter_1', unit='random Gaussian scatter used loc=0, scale=1', format='D', array=scatter_1), fits.Column(name='coolness', unit='0:disturbed, 1:relaxed', format='D', array=coolness)
    ])

    return hdu_cols


hdu_col_array = []
for baseName in baseNames:
    print(baseName)
    path_2_light_cone = os.path.join(test_dir, baseName + '.fits')
    path_2_coordinate_file = os.path.join(
        test_dir, baseName + '_coordinates.fits')
    path_2_galaxy_file = os.path.join(test_dir, baseName + '_galaxy.fits')
    path_2_CLU_file = os.path.join(test_dir, baseName + '_CLU.fits')
    path_2_out_file = os.path.join(tmp_dir, baseName + '_eRo_CLU.fit')

    hdu_col = get_data(
        path_2_light_cone,
        path_2_coordinate_file,
        path_2_galaxy_file,
        path_2_CLU_file)
    tb_hdu = fits.BinTableHDU.from_columns(hdu_col)
    prihdr = fits.Header()
    prihdr['author'] = 'JC'
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tb_hdu])
    if os.path.isfile(path_2_out_file):
        os.system("rm " + path_2_out_file)
    thdulist.writeto(path_2_out_file)
    print('written', path_2_out_file, time.time() - t0)

# concatenates shells into a single fits catalog for the clusters, for
# each simulation
os.chdir(tmp_dir)
c1 = 'ls *_eRo_CLU.fit > fit_list_eRo_CLU.list'
print(c1)
os.system(c1)
#path_2_list = os.path.join(tmp_dir, 'fit_list_eRo_CLU.list')
c2 = stilts_cmd + """ tcat in=@fit_list_eRo_CLU.list ifmt=fits omode=out ofmt=fits out=""" + path_2_CLU_catalog
print(c2)
os.system(c2)
os.system('rm -rf ' + tmp_dir)
