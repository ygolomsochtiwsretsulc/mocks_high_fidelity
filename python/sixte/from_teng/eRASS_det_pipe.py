#!/usr/bin/env python
import sys
import os
import subprocess
import json
import argparse

Version = 10
__doc__ = """
eFEDS point-source detection pipeline
-------------------------

Changes:

Version 10
	* cp eFEDS_det_pipe.py eRASS_det_pipe.py
	* set image size "size=12000 12000"
Version 9
	* change evtool:size from "16000 9000" to "18000 9000"
	* add ersensmap:likemin=MinDetLike
Version n.3
	* Use two bands h1 (2.3-4) and h2 (4-10) instead of the standard 4 bands. I want to make a hard-band detection/selection.
	* 2.3 keV is a turning-point of the effective area curve
	* Calculate ECF using the (PSF- and vignetting-) uncorrected ARF (one_source_in_eFEDS.uncor.arf) and the simput AGN spectral model with NH=1e20 and z=0.
		ECF_2.3_4.0 = 1.313e+11
		ECF_4.0_10. = 3.174e+10
Version 8
	* add a new mode of energy band selection Version n.3
	* add sixte_flag=yes for ermldet, to compensate the sixte-eSASS PSF difference.
	* change erbox:nruns from 3 to 2 in both erbox runs. Georg found that it slightly improves pnt-src/cluster classification. With nruns=3, more point sources will be merged into large sources
	* change ermldet:extmin from 1.5 to 3. (Georg used 2.5 for cluster detection)
	* change emin,emax from 0.2-10 to 0-100 when making EvtImg_full.fits
	* Start to use ermldet:photon_flag=yes. Hopefully the cost of time leads to a better efficiency.
	* Start to use ersensmap with method=FIT, likemin=5, extflag=no
	* change ermldet:cutrad from 20 to 16 in order to make a test.
Version 7
	change ermldet:extlikemin from 10 to 15
Version 6
	change ermldet:extlikemin from 3 to 10
Version 5
	change ermldet:multrad from 15 to 20 (I think spurious extended sources caused by blending are acceptable. They can be recognized later.)
Version n.2
	Use two bands (0.5-2 and 2-10) instead of the standard 4 bands
Version 4
	Use vignetted exposure map in ermldet
	Calculate ECF using the (PSF- and vignetting-) uncorrected ARF (one_source_in_eFEDS.uncor.arf)
	and the simput AGN spectral model with NH=1e20 and z=0.
		ECF_0.2_0.5 = 8.881e+11
		ECF_0.5_1.0 = 1.359e+12
		ECF_1.0_2.0 = 1.010e+12
		ECF_2.0_10. = 9.162e+10
	Start to use eSASSusers_20190520
Version 3
	Add ECF, as measured using the simput AGN spectral model with NH=1e20 and z=0.
	Make an unvignetted exposure map, and use it in ermldet.
	Instead of eSASSusers_190220, use eSASSdevel, in which a few bugs were fixed recently.
Version 2
	change ermldet:cutrad from 15 to 20
"""


def preclear(*arg):
	for files in arg:
		if isinstance(files, str):
			if os.path.isfile(files):
				os.remove(files)
		elif isinstance(files, list):
			for f in files:
				if os.path.isfile(f):
					os.remove(f)
		else:
			sys.exit('ERROR: %s' % files)


Pars2save = ['Version', 'emin_keV', 'emax_keV', 'bandname', 'ecf',
	'band2use', 'infile', 'outdir', 'outprefix', 'cmd_evtool',
	'cmd_expmap', 'cmd_expmap0', 'cmd_ermask', 'cmd_erbox1',
	'cmd_backmap', 'cmd_erbox2', 'cmd_ermldet', 'cmd_sensmap_s', 'cmd_sensmap_h', 'cmd_catprep']

parser = argparse.ArgumentParser(description=__doc__, epilog="")
parser.add_argument('--bands', type=int, default=1, help='index of energy bands setting')
#parser.add_argument('--simulation', action='store_true', help='Simulate??') # obsolete, used to overcome bugs which are now fixed.
parser.add_argument('infile', type=str, help='input event file')
parser.add_argument('outdir', nargs='?', default='', type=str, help='output directory')
args = parser.parse_args()

if args.bands == 1:
	VersApp = ''  # standard 4 band source detection
elif args.bands == 2:
	VersApp = '.2'	# use only two bands
elif args.bands == 3:
	VersApp = '.3'	# use only two bands
else:
	assert False, ('wrong value in --bands', args.bands)

emin_keV = [0.2, 0.5, 1.0, 2.0, 0.5, 2.3, 4.0]
emax_keV = [0.5, 1.0, 2.0, 10., 2.0, 4.0, 10.]
bandname = ["1", "2", "3", "4", "5", "h1", "h2"]
ecf   =    [8.881e+11, 1.359e+12, 1.010e+12, 9.160e+10, 1.178e+12, 1.313e+11, 3.174e+10]
ecf_soft = [3.231e+11, 6.271e+11, 5.440e+11, 1.421e+11, 1.178e+12, 6.297e+10, 3.036e+10]
ecf_hard = [2.083e+11, 4.042e+11, 3.507e+11, 9.160e+10, 7.594e+11, 4.059e+10, 1.956e+10]

band2use = [0, 1, 2, 3]
if args.bands == 2:
	band2use = [4, 3]
elif args.bands == 3:
	band2use = [5, 6]

eFEDSversion = '3'	# Version of eFEDS simulation
outprefix = ""
MinDetLike=8

infile = args.infile
outdir = args.outdir
if outdir=='': outdir = f'./PipeV{Version:d}{VersApp:s}_results'

if not os.path.isfile(infile):
	sys.exit('ERROR: No such file "{}"'.format(infile))

if not os.path.isdir(outdir):
	os.mkdir(outdir)

emin_keV = [emin_keV[i] for i in band2use]
emax_keV = [emax_keV[i] for i in band2use]
bandname = [bandname[i] for i in band2use]
ecf = [ecf[i] for i in band2use]
ecf_soft = [ecf_soft[i] for i in band2use]
ecf_hard = [ecf_hard[i] for i in band2use]

EvtImgFiles = [os.path.join(outdir, f"{outprefix}02{bname}_EvtImg.fits") for bname in bandname]
SrcMapFiles = [os.path.join(outdir, f"{outprefix}02{bname}_SrcMap.fits") for bname in bandname]
ExpMapFiles = [os.path.join(outdir, f"{outprefix}02{bname}_ExpMap.fits") for bname in bandname]
BkgMapFiles = [os.path.join(outdir, f"{outprefix}02{bname}_BkgImg.fits") for bname in bandname]
CheMskFiles = [os.path.join(outdir, f"{outprefix}02{bname}_CheMsk.fits") for bname in bandname]
UnVigExpMap = os.path.join(outdir, f"{outprefix}020_ExpMap.fits")
DetMask = os.path.join(outdir, f"{outprefix}020_DetMsk.fits")
BoxCat1 = os.path.join(outdir, f"{outprefix}020_BoxCa1.fits")
BoxCat2 = os.path.join(outdir, f"{outprefix}020_BoxCa2.fits")
MLCat = os.path.join(outdir, f"{outprefix}020_MLCat.fits")
SrcCat = os.path.join(outdir, f"{outprefix}020_SrcCat.fits")
SkyCovSoft = os.path.join(outdir, f"{outprefix}SkyCov_Soft_{MinDetLike:g}.fits")
SenMapSoft = os.path.join(outdir, f"{outprefix}SenMap_Soft_{MinDetLike:g}.fits")
SkyCovHard = os.path.join(outdir, f"{outprefix}SkyCov_Hard_{MinDetLike:g}.fits")
SenMapHard = os.path.join(outdir, f"{outprefix}SenMap_Hard_{MinDetLike:g}.fits")

cmd_evtool = []
if 1:
	cmd_evtool.append(["evtool",
	f"eventfiles={infile}",
	f"outfile={os.path.join(outdir,'EvtImg_full.fits')}",
	f"emin=0",
	f"emax=100",
	f"image=yes",
	f"rebin=80",
	f"size=12000 12000",
	f"pattern=15",
	f"center_position=0 0"
	])
for i in range(len(bandname)):
	cmd_evtool.append(["evtool",
	f"eventfiles={infile}",
	f"outfile={EvtImgFiles[i]}",
	f"emin={emin_keV[i]:f}",
	f"emax={emax_keV[i]:f}",
	f"image=yes",
	f"rebin=80",
	f"size=12000 12000",
	f"pattern=15",
	f"center_position=0 0"
	])

cmd_expmap = ["expmap",
	f"inputdatasets={infile}",
	f"emin=" + " ".join(f"{x:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x:g}" for x in emax_keV),
	f"templateimage={EvtImgFiles[0]}",
	f"withvignetting=yes",
	f"withmergedmaps=yes",
	f"withfilebadpix=no",
	f"withcalbadpix=no",
	f"mergedmaps={' '.join(ExpMapFiles)}"
	]

cmd_expmap0 = ["expmap",
	f"inputdatasets={infile}",
	f"emin=0.5",
	f"emax=2.0",
	f"templateimage={EvtImgFiles[0]}",
	f"withvignetting=no",
	f"withmergedmaps=yes",
	f"withfilebadpix=no",
	f"withcalbadpix=no",
	f"mergedmaps={UnVigExpMap}"
	]

cmd_ermask = ["ermask",
	f"expimage={ExpMapFiles[0]}",
	f"detmask={DetMask}",
	f"threshold1=0.2",
	f"threshold2=1.0",
	f"regionfile_flag=no"
	]

cmd_erbox1 = ["erbox",
	f"images={' '.join(EvtImgFiles)}",
	f"boxlist={BoxCat1}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"detmasks={DetMask}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"hrdef=",
	f"ecf=" + " ".join(f"{x:g}" for x in ecf),
	f"nruns=2",
	f"likemin=6.0",
	f"boxsize=4",
	f"compress_flag=N",
	f"bkgima_flag=N",
	f"expima_flag=Y",
	f"detmask_flag=Y"
	]

cmd_backmap = []
for i in range(len(band2use)):
	cmd_backmap.append(["erbackmap",
	f"image={EvtImgFiles[i]}",
	f"expimage={ExpMapFiles[i]}",
	f"boxlist={BoxCat1}",
	f"detmask={DetMask}",
	f"cheesemask={CheMskFiles[i]}",
	f"bkgimage={BkgMapFiles[i]}",
	f"idband={i:d}",
	f"scut=0.0001",
	f"mlmin=0",
	f"maxcut=0.5",
	f"fitmethod=smooth",
	f"nsplinenodes=36",
	f"degree=2",
	f"smoothflag=yes",
	f"smoothval=15.",
	f"snr=40.0",
	f"excesssigma=1000.",
	f"nfitrun=3",
	f"cheesemaskflag=Y"
	])

cmd_erbox2 = ["erbox",
	f"images={' '.join(EvtImgFiles)}",
	f"boxlist={BoxCat2}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"detmasks={DetMask}",
	f"bkgimages={' '.join(BkgMapFiles)}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"hrdef=",
	f"ecf=" + " ".join(f"{x:g}" for x in ecf),
	f"nruns=2",
	f"likemin=4.",
	f"boxsize=4",
	f"compress_flag=N",
	f"bkgima_flag=Y",
	f"expima_flag=Y",
	f"detmask_flag=Y"
	]

cmd_ermldet = ["ermldet",
	f"mllist={MLCat}",
	f"boxlist={BoxCat2}",
	f"images={' '.join(EvtImgFiles)}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"detmasks={DetMask}",
	f"bkgimages={' '.join(BkgMapFiles)}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"hrdef=",
	f"ecf=" + " ".join(f"{x:g}" for x in ecf),
	f"likemin=5.",
	f"extlikemin=15",
	f"compress_flag=N",
	f"cutrad=16",
	f"multrad=20",
	f"extmin=2.5",
	f"extmax=30.0",
	f"bkgima_flag=Y",
	f"expima_flag=Y",
	f"detmask_flag=Y",
	f"extentmodel=beta",
	f"thres_flag=N",
	f"thres_col=like",
	f"thres_val=30.",
	f"nmaxfit=3",
	f"nmulsou=2",
	f"fitext_flag=yes",
	f"srcima_flag=yes",
	f"srcimages={' '.join(SrcMapFiles)}",
	f"shapelet_flag=yes",
	f"photon_flag=no",
	f"sixte_flag=yes"
	]

cmd_catprep = ["catprep", f"infile={MLCat}", f"outfile={SrcCat}"]

cmd_sensmap_s = ["ersensmap",
	f"sensimage={SenMapSoft}",
	f"area_table={SkyCovSoft}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"detmasks={DetMask}",
	f"bkgimages={' '.join(BkgMapFiles)}",
	f"srcimages={' '.join(SrcMapFiles)}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"ecf=" + " ".join(f"{x:g}" for x in ecf_soft),
	f"method=FIT",
	f"likemin={MinDetLike:g}",
	f"extflag=no",
	f"extentmodel=beta",
	f"extlikemin=15",
	f"detmask_flag=yes",
	f"shapelet_flag=no",
	f"photon_flag=no",
	f"area_flag=yes"
	]

cmd_sensmap_h = ["ersensmap",
	f"sensimage={SenMapHard}",
	f"area_table={SkyCovHard}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"detmasks={DetMask}",
	f"bkgimages={' '.join(BkgMapFiles)}",
	f"srcimages={' '.join(SrcMapFiles)}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"ecf=" + " ".join(f"{x:g}" for x in ecf_hard),
	f"method=FIT",
	f"likemin={MinDetLike:g}",
	f"extflag=no",
	f"extentmodel=beta",
	f"extlikemin=15",
	f"detmask_flag=yes",
	f"shapelet_flag=no",
	f"photon_flag=no",
	f"area_flag=yes"
	]

################################################################################
if 1:
	localvars = locals()
	logfile = os.path.join(outdir, f'eFEDS_V{Version:03d}{VersApp}.par')
	#if os.path.isfile(logfile): sys.exit(f'ERROR: {logfile} already exists, please delete it')
	Parsdict = dict((k, localvars[k]) for k in Pars2save)
	Parsdict['History'] = __doc__.split('\n')
	with open(logfile, 'w') as flog:
		json.dump(Parsdict, flog, indent=4)

if 1:
	preclear(EvtImgFiles)
	for cmd in cmd_evtool:
		subprocess.run(cmd, check=True)

	preclear(ExpMapFiles)
	subprocess.run(cmd_expmap, check=True)

	preclear(UnVigExpMap)
	subprocess.run(cmd_expmap0, check=True)

	preclear(DetMask)
	subprocess.run(cmd_ermask, check=True)

	preclear(BoxCat1)
	subprocess.run(cmd_erbox1, check=True)

	preclear(BkgMapFiles, CheMskFiles)
	for cmd in cmd_backmap:
		subprocess.run(cmd, check=True)

	preclear(BoxCat2)
	subprocess.run(cmd_erbox2, check=True)

	preclear(SrcMapFiles, MLCat)
	subprocess.run(cmd_ermldet, check=True)

	preclear(SrcCat)
	subprocess.run(cmd_catprep, check=True)

	preclear(SenMapSoft, SkyCovSoft)
	subprocess.run(cmd_sensmap_s, check=True)
