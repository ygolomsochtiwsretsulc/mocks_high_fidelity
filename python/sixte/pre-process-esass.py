"""
Preprocess the simulated event files to make compatible with eSASS.
A. Malyali, 2019. amalyali@mpe.mpg.de

source /home/erosita/sw/sass-setup.sh eSASSusers_190520


import os, sys

tile_id = '355'
folder = "/data40s/erosim/eRASS/extendedSourceCat-v0.1/" + tile_id
os.mkdir(folder)

cmd = "/data40s/erosim/eRASS/extendedSourceCat-v0.1/" + tile_id
os.chdir(cmd)

list_cmd = "ln -s /data40s/erosim/eRASS/sixte/" + tile_id +"/* ."
os.system(list_cmd)

List of commands run 
7 times ero_calevents
7 times fparkey
1 time evtool
1 time radec2xy

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd1_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd1_evt.fits CCDNr=1 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd1_evt.fits[1] keyword=PAT_SEL value=15 add=yes

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd2_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd2_evt.fits CCDNr=2 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd2_evt.fits[1] keyword=PAT_SEL value=15 add=yes

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd3_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd3_evt.fits CCDNr=3 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd3_evt.fits[1] keyword=PAT_SEL value=15 add=yes

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd4_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd4_evt.fits CCDNr=4 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd4_evt.fits[1] keyword=PAT_SEL value=15 add=yes

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd5_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd5_evt.fits CCDNr=5 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd5_evt.fits[1] keyword=PAT_SEL value=15 add=yes

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd6_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd6_evt.fits CCDNr=6 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd6_evt.fits[1] keyword=PAT_SEL value=15 add=yes

ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd7_evt.fits eroEvtFile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd7_evt.fits CCDNr=7 RA=0.0 Dec=-4.780191847199159
fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd7_evt.fits[1] keyword=PAT_SEL value=15 add=yes

evtool eventfiles="/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd1_evt.fits /data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd2_evt.fits /data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd3_evt.fits /data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd4_evt.fits /data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd5_evt.fits /data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd6_evt.fits /data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd7_evt.fits" outfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/merged_erass.fits pattern=15 clobber=yes

radec2xy file=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/merged_erass.fits ra0=0.0 dec0=-4.780191847199159


ero_calevents Projection=AIT Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits clobber=yes EvtFile=/data40s/erosim/eRASS/eRASS8/400/erass_ccd7_evt.fits eroEvtFile=/data40s/erosim/eRASS/test.fits CCDNr=7 RA=0.0 Dec=-4.780191847199159 RA_PNT=0.0 DEC_PNT=-4.780191847199159 

fparkey fitsfile=/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/400/cal_erass_ccd7_evt.fits[1]

"""
import subprocess
import os
import sys
from astropy_healpix import healpy
import numpy as n
pix_ids = n.arange(healpy.nside2npix(8) )
ra_cen_s, dec_cen_s = healpy.pix2ang(8, pix_ids, nest=False, lonlat=True)

ccds = [1, 2, 3, 4, 5, 6, 7]


class PrepareForEsass:
    def __init__(self, ra_cen, dec_cen, data_dir, out_dir):
        """
        :param ra_cen: Right ascension of skyfield center position
        :param dec_cen: Declination of skyfield center position
        directory with inputs, directory with outputs
        """
        self._data_dir = data_dir
        self._out_dir = out_dir
        self._ra_cen = ra_cen
        self._dec_cen = dec_cen
        try:
            os.makedirs(self._out_dir)
        except OSError as e:
            print('already created')
    def cal_events(self):
        """
        Prepare event lists from SIXTE simulation for the eSASS pipeline.
        Calibrate event lists from SIXTE (mainly to ensure we have the correct extensions for the FITS files).
        """
        for ii in ccds:
            uncal_file = '%s/erass_ccd%s_evt.fits' % (self._data_dir, ii)
            cal_file = '%s/cal_erass_ccd%s_evt.fits' % (self._out_dir, ii)

            cmd = ["ero_calevents",
                   "Projection=AIT",
                   "Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits",
                   "clobber=yes",
                   "EvtFile=%s" % uncal_file,
                   "eroEvtFile=%s" % cal_file,
                   "CCDNr=%i" % ii,
                   "RA=%s" % self._ra_cen,
                   "Dec=%s" % self._dec_cen  # center of skyfield
                   ]
            print(" ".join(cmd))
            subprocess.check_call(cmd)

            cmd = ["fparkey",
                   "fitsfile=%s[1]" % cal_file,
                   "keyword=PAT_SEL",
                   "value=15",
                   "add=yes"
                   ]
            print(" ".join(cmd))
            subprocess.check_call(cmd)

    def merge_cal_events_across_ccds(self):
        """
        1. Merge all calibrated event files from all CCDs.
        2. Centre ra dec 2 xy
        """
        unmerged_cal_evt_files = []
        for ii in ccds:
            unmerged_cal_evt_files.append(
                '%s/cal_erass_ccd%s_evt.fits' %
                (self._out_dir, ii))

        merged_cal_evt_file = "%s/merged_erass.fits" % self._out_dir

        cmd = ["evtool", "eventfiles=%s" %
               ' '.join([str(x) for x in unmerged_cal_evt_files]), "outfile=%s" %
               merged_cal_evt_file, "pattern=15", "clobber=yes"]
        print(" ".join(cmd))
        subprocess.check_call(cmd)

        cmd = ["radec2xy",
               "file=%s" % merged_cal_evt_file,
               "ra0=%s" % self._ra_cen,
               "dec0=%s" % self._dec_cen
               ]
        print(" ".join(cmd))
        subprocess.check_call(cmd)

    def prep_and_merge(self):
        self.cal_events()
        self.merge_cal_events_across_ccds()


if __name__ == '__main__':
    # Define parameters for simulation
    tile_id = sys.argv[1]  # '355'
    ra_cen = ra_cen_s[int(tile_id)]
    dec_cen = dec_cen_s[int(tile_id)]
    data_dir = "/data40s/erosim/eRASS/eRASS8/" + tile_id
    out_dir = "/data40s/erosim/eRASS/extendedSourceCat-eRASS8-v0.1/" + tile_id
    # Preprocess event files...
    PrepareForEsass(ra_cen, dec_cen, data_dir, out_dir).prep_and_merge()
