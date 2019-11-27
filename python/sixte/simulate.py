"""
Simulate eROSITA event data from sky model using SIXTE.
A. Malyali, 2019. amalyali@mpe.mpg.de

modification by J.C:
 - more flexible with directories and multiple simput files
 - merge of the 3 gti files obtained using 'mgtime'

Full sky run :
On ds52 only :
nohup sh xaa > xaa.log &
nohup sh xab > xab.log &
nohup sh xac > xac.log &
nohup sh xad > xad.log &
nohup sh xae > xae.log &
nohup sh xaf > xaf.log &
nohup sh xag > xag.log &
nohup sh xah > xah.log &

ero_vis GTIfile=/data40s/erosim/eRASS/sixte/000/erass.gti Simput=/data17s/darksim/MD/MD_1.0Gpc/cat_AGN_SIMPUT/SIMPUT_000000_1024.fit Exposure=31536000.000000 Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits TSTART=0.000000 dt=1.0 visibility_range=1.0 clobber=yes

source /home/erosita/sw/sass-setup.sh eSASSdevel 

ls /data40s/erosim/eRASS/sixte/*/erass_ccd?_evt.fits > list_ccd1.lis
stilts tcat in=@list_ccd1.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/simulated_photons.fits ofmt=fits

cp /data40s/erosim/eRASS/simulated_photons.fits wwwDir/erosita_stuff/

cd data/eRoMok
wget http://www.mpe.mpg.de/~comparat/erosita_stuff/simulated_photons.fits

"""
import subprocess
import os
import errno
import sys
from astropy_healpix import healpy
import numpy as n
pix_ids = n.arange(healpy.nside2npix(8) )
ra_cen_s, dec_cen_s = healpy.pix2ang(8, pix_ids, nest=False, lonlat=True)


class Simulator:
    """
    SIXTE simulator for eROSITA observations.
    1. Compute GTI file for given simput
    2. Simulate eROSITA observations of simput, using GTI to speed things up.
    """

    def __init__(self, with_bkg_par, t_start, exposure, seed, simput, simput2, simput3, data_dir, ra_cen, dec_cen):
        # def __init__(self, with_bkg_par, t_start, exposure, seed, simput,
        # simput2, simput3):
        """
        :param with_bkg_par: Simulate with particle background.
        :param t_start: Start time of simulation. Input units of [s]
        :param exposure: Length of time to simulate for after t_start
        :param seed: Seed for random number generator.
        :param simput: Simput file (ie. the sky model)
        """
        self._with_bkg_par = bool(with_bkg_par)
        self._t_start = float(t_start)  # secs
        self._exposure = float(exposure)
        self._seed = int(seed)
        self._simput = simput
        self._simput2 = simput2
        self._simput3 = simput3
        self._data_dir = data_dir
        self._ra_cen = ra_cen
        self._dec_cen = dec_cen

    def make_event_directory(self):
        """
        Check for whether directory exists and create if not.
        """
        try:
            os.makedirs(self._data_dir)
        except OSError as e:
            print('already exists')

    def compute_gti(self, simput_ero_vis, ):
        """
        Compute the GTI (good time interval) file given the simput file.
        Use this as input to SIXTE call for reducing computational expense.
        """
        gti_file = os.path.join(
            self._data_dir,
            "erass_" +
            os.path.basename(simput_ero_vis)[
                :-
                4] +
            ".gti")

        cmd = ["ero_vis",
               "GTIfile=%s" % gti_file,
               "Simput=%s" % simput_ero_vis,
               "Exposure=%f" % self._exposure,
               "Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits",
               "TSTART=%f" % self._t_start,
               #"RA=%s" % self._ra_cen,
               #"Dec=%s" % self._dec_cen,
               "dt=1.0",
               "visibility_range=1.0",
               "clobber=yes"
               ]

        print(" ".join(cmd))
        subprocess.check_call(cmd)

    def run_sixte(self):
        """
        Launch erosim from python.
        """
        prefix = "%s/erass_" % self._data_dir

        cmd = ["erosim",
               "Simput=%s" % self._simput,
               "Simput2=%s" % self._simput2+ "[IMGSCAL>0]",
               "Simput3=" + self._simput3,
               "Prefix=%s" % prefix,
               "Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits",
               "RA=%s" % self._ra_cen,
               "Dec=%s" % self._dec_cen,
               "GTIFile=%s/erass.gti" % self._data_dir,
               "TSTART=%s" % self._t_start,
               "Exposure=%s" % self._exposure,
               "MJDREF=51543.875",
               "dt=1.0",
               "Seed=%s" % self._seed,
               "clobber=yes",
               "chatter=3"
               ]

        if self._with_bkg_par is True:
            cmd.append("Background=yes")
        else:
            cmd.append("Background=no")

        print(" ".join(cmd))
        subprocess.check_call(cmd)

    def run_all(self):
        """
        Run SIXTE simulation of eRASS 1
        """
        print('make event directory')
        self.make_event_directory()
        print('compute gti with ero_vis')
        self.compute_gti(self._simput)
        self.compute_gti(self._simput2)
        self.compute_gti(self._simput3)

        print('merges gti files')
        path_to_gti_list = os.path.join(self._data_dir, "gti.list")
        #print(path_to_gti_list)
        path_to_gti = os.path.join(self._data_dir, "erass.gti")
        print(path_to_gti)
        command_list = "ls " + self._data_dir + "/*.gti > " + path_to_gti_list
        print(command_list)
        os.system(command_list)
        # ls /data40s/erosim/eRASS/sixte/000/erass_SIMPUT_000000*.gti >
        # gti.list
        command_merge = "mgtime ingtis=@" + path_to_gti_list + \
            " outgti=" + path_to_gti + " merge=OR"
        print(command_merge)
        os.system(command_merge)

        print('run sixte with erosim')
        self.run_sixte()


if __name__ == '__main__':
    # Define parameters for simulation
    bkg = 1
    t_start = 0.0
    exposure = 31536000 * 4 # = 4 years  # 31536000 = 1year # 15750000 = 1/2 year
    seed = 42
    tile_id = sys.argv[1]  # '355'
    ra_cen = ra_cen_s[int(tile_id)]
    dec_cen = dec_cen_s[int(tile_id)]
    data_dir = "/data40s/erosim/eRASS/eRASS8/" + tile_id
    simput_file_agn = '/data17s/darksim/MD/MD_1.0Gpc/cat_AGN_SIMPUT/SIMPUT_000' + tile_id + '_1024.fit'
    simput_file_cluster = '/data17s/darksim/MD/MD_1.0Gpc/cat_CLU_SIMPUT/SIMPUT_000' + tile_id + '.fit'
    simput_file_BKG = '/data40s/erosim/Simput/Fullsky_bkg/bkg_simput_' + str(int(tile_id)) + '.fits'
    # Launch...
    # 3 files
    #Simulator(bkg, t_start, exposure, seed, simput_file_1, simput_file, simput_file_2).run_all()
    # 2 files
    Simulator(
        bkg,
        t_start,
        exposure,
        seed,
        simput_file_agn,
        simput_file_cluster,
        simput_file_BKG, 
        data_dir,
        ra_cen,
        dec_cen).run_all()

# /data17s/darksim/MD/MD_1.0Gpc/cat_AGN_SIMPUT/SIMPUT_000355_1024.fit
