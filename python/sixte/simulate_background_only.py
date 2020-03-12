"""
Simulate eROSITA background events with the model from T. Boller, adapted by T. Liu

ero_vis GTIfile=/data40s/erosim/eRASS/background/009/erass.gti Simput=/data40s/erosim/Simput/Fullsky_bkg/bkg_simput_9.fits Exposure=31536000.000000 Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits TSTART=0.000000 dt=10.0 visibility_range=1.0 clobber=yes

ero_vis GTIfile=/data40s/erosim/eRASS/background/009/erass.gti Simput=/data17s/darksim/MD/MD_1.0Gpc/cat_CLU_SIMPUT/SIMPUT_000009.fit Exposure=31536000.000000 Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits TSTART=0.000000 dt=10.0 visibility_range=1.0 clobber=yes

erosim Simput=/data40s/erosim/Simput/Fullsky_bkg/bkg_simput_9.fits Prefix=/data40s/erosim/eRASS/background/009/erass_ Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits GTIFile=/data40s/erosim/eRASS/background/009/erass.gti TSTART=0.0 Exposure=31536000.0 MJDREF=51543.875 dt=10.0 Seed=42 clobber=yes chatter=3 Background=no

erosim Simput=/data17s/darksim/MD/MD_1.0Gpc/cat_CLU_SIMPUT/SIMPUT_000009.fit Prefix=/data40s/erosim/eRASS/background/009/erass_ Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits GTIFile=/data40s/erosim/eRASS/background/009/erass.gti TSTART=0.0 Exposure=31536000.0 MJDREF=51543.875 dt=10.0 Seed=42 clobber=yes chatter=3 Background=no

erosim Simput=bkg_simput_9.fits Prefix=erass_ Attitude=eRASS_Pc87M55_3dobi_att_remeis.fits GTIFile=erass.gti TSTART=0.0 Exposure=31536000.0 MJDREF=51543.875 dt=10.0 Seed=42 clobber=yes chatter=3 Background=no

erosim Simput=SIMPUT_000009.fit Prefix=erass_ Attitude=eRASS_Pc87M55_3dobi_att_remeis.fits GTIFile=erass.gti TSTART=0.0 Exposure=31536000.0 MJDREF=51543.875 dt=10.0 Seed=42 clobber=yes chatter=3 Background=no


#cd /data40s/erosim/eRASS/

"""
import subprocess
import os
import errno
import sys


class Simulator:
    """
    SIXTE simulator for eROSITA observations.
    1. Compute GTI file for given simput
    2. Simulate eROSITA observations of simput, using GTI to speed things up.
    """

    def __init__(
            self,
            with_bkg_par,
            t_start,
            exposure,
            seed,
            simput,
            data_dir, tile_id):
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
        self._data_dir = data_dir
        self._tile_id = tile_id

    def make_event_directory(self):
        """
        Check for whether directory exists and create if not.
        """
        try:
            os.makedirs(self._data_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def compute_gti(self, simput_ero_vis, ):
        """
        Compute the GTI (good time interval) file given the simput file.
        Use this as input to SIXTE call for reducing computational expense.
        """
        gti_file = os.path.join(
            self._data_dir, "erass.gti")

        cmd = ["ero_vis",
               "GTIfile=%s" % gti_file,
               "Simput=%s" % simput_ero_vis,
               "Exposure=%f" % self._exposure,
               "Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits",
               "TSTART=%f" % self._t_start,
               # "RA=0.0",
               # "Dec=0.0",
               "dt=10.0",
               "visibility_range=1.0",
               "clobber=yes"
               ]

        #subprocess.check_call(cmd)
        #print(cmd)
        cmd = "cp /data40s/erosim/eRASS/sixte/"+self._tile_id.zfill(3)+"/erass.gti /data40s/erosim/eRASS/background/"+self._tile_id.zfill(3)+"/erass.gti"
        print(cmd)
        os.system(cmd)


    def run_sixte(self):
        """
        Launch erosim from python.
        """
        prefix = "%s/erass_" % self._data_dir

        cmd = ["erosim",
               "Simput=%s" % self._simput,
               "Prefix=%s" % prefix,
               "Attitude=/data40s/erosim/eRASS/eRASS_Pc87M55_3dobi_att_remeis.fits",
               # "RA=0.0",
               # "Dec=0.0",
               "GTIFile=%s/erass.gti" % self._data_dir,
               "TSTART=%s" % self._t_start,
               "Exposure=%s" % self._exposure,
               "MJDREF=51543.875",
               "dt=10.0",
               "Seed=%s" % self._seed,
               "clobber=yes",
               "chatter=3"
               ]

        if self._with_bkg_par is True:
            cmd.append("Background=yes")
        else:
            cmd.append("Background=no")

        subprocess.check_call(cmd)
        print(cmd)


    #def cmd_postprocess(self):
        #ccd1_file = "/data40s/erosim/eRASS/background/"+self._tile_id"/erass_ccd1_evt.fits"
        #out_ccd1_file = "/data40s/erosim/eRASS/background/"+self._tile_id"/softbkg_ccd1_evt.fits"
        #command = "stilts tpipe in="+ccd1_file+""" ifmt=fits icmd=explodeall cmd='addcol HPX " healpixRingIndex( 3, RA, -DEC )"' \ cmd='select "HPX=="""+self._tile_id+""""' omode=out ofmt=fits out="""+out_ccd1_file
        #return command
    
    def run_all(self):
        """
        Run SIXTE simulation of eRASS 1
        """
        print('make event directory')
        self.make_event_directory()
        print('copy gti with ero_vis')
        self.compute_gti(self._simput)
        print('run sixte with erosim')
        self.run_sixte()
        #cmd = self.cmd_postprocess()
        #os.system(cmd)


if __name__ == '__main__':
    # Define parameters for simulation
    bkg = 0
    t_start = 0.0
    exposure = 31536000  # =1year # 15750000 # = 1/2 year
    seed = 42
    tile_id = sys.argv[1]  # '355'
    data_dir = "/data40s/erosim/eRASS/background/" + tile_id.zfill(3)
    simput_file_BKG = '/data40s/erosim/Simput/Fullsky_bkg/bkg_simput_' + tile_id + '.fits'

    Simulator(
        bkg,
        t_start,
        exposure,
        seed,
        simput_file_BKG,
        data_dir, tile_id).run_all()
