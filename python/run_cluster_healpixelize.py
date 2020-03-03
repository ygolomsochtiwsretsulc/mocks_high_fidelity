"""
export GIT_AGN_MOCK='/home/comparat/software/lss_mock_dev/'
export GIT_XPAINT='/home/comparat/software/XSB_Painting'


"""
import glob
import sys
from astropy_healpix import healpy
import os
from astropy.table import Table
import astropy.io.fits as fits
import numpy as n
import time
t0 = time.time()

env = sys.argv[1]  # 'MD10'

test_dir = os.path.join(os.environ[env], 'fits')
catalogue_CLU_dir = os.path.join(os.environ[env], 'cat_eRO_CLU')
catalogue_RS_dir = os.path.join(os.environ[env], 'cat_eRO_CLU_RS')

if os.path.isdir(catalogue_CLU_dir) == False:
    os.system('mkdir -p ' + catalogue_CLU_dir)
if os.path.isdir(catalogue_RS_dir) == False:
    os.system('mkdir -p ' + catalogue_RS_dir)

path_2_eRO_CLU = os.path.join(test_dir, env + '_eRO_CLU.fit')
path_2_eRO_CLU_RS = os.path.join(test_dir, env + '_eRO_CLU_SAT_RS.fit')

lss_git_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python')


def split_CLU():
    # eRO_CLU
    hd = fits.open(path_2_eRO_CLU)

    HEALPIX_32 = healpy.ang2pix(
        8,
        hd[1].data['dec'] *
        n.pi /
        180. +
        n.pi /
        2.,
        hd[1].data['ra'] *
        n.pi /
        180.)

    for HEALPIX_32_id in n.arange(healpy.nside2npix(8)):
        sf = (HEALPIX_32 == HEALPIX_32_id)
        path_2_agn_summary_file = os.path.join(
            catalogue_CLU_dir, str(HEALPIX_32_id).zfill(6) + '.fit')
        if os.path.isfile(path_2_agn_summary_file):
            os.system('rm ' + path_2_agn_summary_file)
        t = Table(hd[1].data.T[sf])
        t.write(path_2_agn_summary_file, format='fits')

    hd.close()


def split_CLU_RS():
    # eRO_CLU_RS
    hd = fits.open(path_2_eRO_CLU_RS)

    HEALPIX_32 = healpy.ang2pix(
        8,
        hd[1].data['dec'] *
        n.pi /
        180. +
        n.pi /
        2.,
        hd[1].data['ra'] *
        n.pi /
        180.)

    for HEALPIX_32_id in n.arange(healpy.nside2npix(8)):
        sf = (HEALPIX_32 == HEALPIX_32_id)
        path_2_agn_summary_file = os.path.join(
            catalogue_RS_dir, str(HEALPIX_32_id).zfill(6) + '.fit')
        if os.path.isfile(path_2_agn_summary_file):
            os.system('rm ' + path_2_agn_summary_file)
        t = Table(hd[1].data.T[sf])
        t.write(path_2_agn_summary_file, format='fits')

    hd.close()


split_CLU()

split_CLU_RS()
