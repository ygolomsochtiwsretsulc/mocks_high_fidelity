#!/usr/bin/env python
import sys, os
import subprocess
import json
import argparse
from eRO_settings import Bands
from astropy.io import fits

Version = 1
__doc__ = """
prepare images and exposure maps

Written by T. Liu (MPE)

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
	'cmd_expmap', 'cmd_expmap0', 'cmd_ermask']

parser = argparse.ArgumentParser(description=__doc__, epilog="")
parser.add_argument('infile', type=str, help='input event file')
parser.add_argument('outdir', type=str, help='output directory')
parser.add_argument('-b', dest='band', type=int, default=-1, help='index of energy band')
parser.add_argument('-C160','--C160',dest='C160', action='store_true')
parser.add_argument('-C180','--C180',dest='C180', action='store_true')
parser.add_argument('-SEP','--SEP',dest='SEP', action='store_true')
parser.add_argument('-hp8','--hp8',dest='hp8', action='store_true')
args = parser.parse_args()

if args.C180:
	os.environ["CALDB"]='/data11s/lewton/eFEDS/caldb_20200228_C180'
	print(f'CALDB={os.environ["CALDB"]}')
elif args.C160:
	os.environ["CALDB"]='/data11s/lewton/eFEDS/caldb_C160'
	print(f'CALDB={os.environ["CALDB"]}')

ImgSize='18000 9000'
#ImgSize='9000 9000'
ImgRebin=80
PATTERN=15
#ImgRebin=800
#PATTERN=3

if args.SEP: ImgSize='3600 3600'
if args.hp8:
	ImgSize='12000 12000'
	ImgRebin=80

infile = args.infile
outdir = args.outdir
if not os.path.isdir(outdir): os.mkdir(outdir)
if not os.path.isfile(infile): sys.exit('ERROR: No such file "{}"'.format(infile))

if fits.getval(infile,'CREATOR',ext=1)=='SIXTE': Sixte=True
else: Sixte=False

if Sixte: evflag=0
else: evflag="0xc00fff30"

healpix_id = outdir.split('/')[-2]
outprefix = healpix_id + "_" # ""
print(outprefix)

if args.band==-1:
	band2use = [9,11,12] #,13,4,5]
else:
	assert args.band in Bands.keys()
	band2use=[args.band]

emin_keV = [Bands[i][0] for i in band2use]
emax_keV = [Bands[i][1] for i in band2use]
bandname = [str(i) for i in band2use]
ecf = [Bands[i][2] for i in band2use]


EvtImgFiles = [os.path.join(outdir, f"{outprefix}02{bname}_EvtImg.fits") for bname in bandname]
ExpMapFiles = [os.path.join(outdir, f"{outprefix}02{bname}_ExpMap.fits") for bname in bandname]
CheMskFiles = [os.path.join(outdir, f"{outprefix}02{bname}_CheMsk.fits") for bname in bandname]
UnVigExpMap = [os.path.join(outdir, f"{outprefix}02{bname}_UnvExpMap.fits") for bname in bandname]
DetMask     = [os.path.join(outdir, f"{outprefix}02{bname}_DetMsk.fits") for bname in bandname]

# commands for single band run 
# for multiple bands, add a loop and change [0] for [i]

#cmd_evtool = []
#for i in range(len(bandname)):
cmd_evtool = ["evtool",
	f"eventfiles={infile}", f"outfile={EvtImgFiles[0]}",
	f"emin={emin_keV[0]:f}", f"emax={emax_keV[0]:f}",
	f"events=yes", f"image=yes", f"rebin={ImgRebin:d}",
	f"size={ImgSize}", f"pattern={PATTERN}", f"center_position=0 0",
	f"flag={evflag}"#, "GTI='617943605 649479605' "
	]

cmd_expmap = ["expmap",
	f"inputdatasets={infile}",
	f"emin=" + " ".join(f"{x:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x:g}" for x in emax_keV),
	f"templateimage={EvtImgFiles[0]}",
	f"withdetmaps=yes",
	f"withvignetting=yes",
	f"withsinglemaps=no",
	f"singlemaps=ccd1.fits ccd2.fits ccd3.fits ccd4.fits ccd5.fits ccd6.fits ccd7.fits",
	f"withmergedmaps=yes",
	f"withfilebadpix=yes",
	f"withcalbadpix=yes",
	f"withweights=yes",
	f"mergedmaps={ExpMapFiles[0]}"
	]

cmd_expmap_unvign = ["expmap",
	f"inputdatasets={infile}",
	f"emin=" + " ".join(f"{x:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x:g}" for x in emax_keV),
	f"templateimage={EvtImgFiles[0]}",
	f"withdetmaps=yes",
	f"withvignetting=no",
	f"withmergedmaps=yes",
	f"withfilebadpix=yes",
	f"withcalbadpix=yes",
	f"withweights=yes",
	f"mergedmaps={UnVigExpMap[0]}"
	]

cmd_ermask = ["ermask",
	f"expimage={ExpMapFiles[0]}",
	f"detmask={DetMask[0]}",
	f"threshold1=0.01",
	f"threshold2=1.0",
	f"regionfile_flag=no"
	]

################################################################################
preclear(EvtImgFiles)
#for cmd in cmd_evtool:
if not os.path.isfile(EvtImgFiles[0]):
	print(cmd_evtool)
	subprocess.run(cmd_evtool, check=True)

preclear(ExpMapFiles)
if not os.path.isfile(ExpMapFiles[0]):
	print(cmd_expmap)
	subprocess.run(cmd_expmap, check=True)

preclear(UnVigExpMap)
if not os.path.isfile(UnVigExpMap[0]):
	print(cmd_expmap_unvign)
	subprocess.run(cmd_expmap_unvign, check=True)

preclear(DetMask)
if not os.path.isfile(DetMask[0]):
	print(cmd_ermask)
	subprocess.run(cmd_ermask, check=True)

os.system("chgrp erosim "+outdir+"/*")