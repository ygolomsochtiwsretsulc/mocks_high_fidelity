#!/usr/bin/env python
import sys
import os
import subprocess
import json
import argparse
from eRO_settings import Bands,Bandgroups

Version = 9

__doc__ = """
point-source detection pipeline
From V9, change from extmin=5 extmax=40 to extmin=2 extmax=15
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

def runcmd(cmd,logfile):
	proc=subprocess.run(cmd,check=False,shell=True,text=True,stdout=None,stderr=subprocess.PIPE)
	if proc.returncode!=0:
		print(proc.stderr)
		with open(logfile, 'a') as flog: flog.write(proc.stderr)

parser = argparse.ArgumentParser(description=__doc__, epilog="")
parser.add_argument('inputdir', type=str, help='input')
parser.add_argument('outdir', nargs='?', default='', type=str, help='output directory')
parser.add_argument('-b', dest='VerBand', type=str, default='default', help='energy bands setting')
parser.add_argument('-V', dest='Version', type=int, default='3', help='Version')
#parser.add_argument('--simulation', action='store_true', help='Simulate??') # obsolete, used to overcome bugs which are now fixed.
parser.add_argument('-C160','--C160',dest='C160', action='store_true')
parser.add_argument('-C180','--C180',dest='C180', action='store_true')
args = parser.parse_args()

if args.C180:
	os.environ["CALDB"]='/data11s/lewton/eFEDS/caldb_20200228_C180'
	print(f'CALDB={os.environ["CALDB"]}')
elif args.C160:
	os.environ["CALDB"]='/data11s/lewton/eFEDS/caldb_C160'
	print(f'CALDB={os.environ["CALDB"]}')

VerBand = args.VerBand
Version = args.Version
#assert Version in [3,4,5]
#assert Version in [6,7,8]
#assert Version in [9,10,11,12]
band2use = Bandgroups[args.VerBand]
if VerBand=='default': VerBand=''

MinDetLike=10
ERML_MULTRAD=20
Sixte=False
if Version==2 or Version==3:
	ERML_CUTRAD=16
	ERML_likemin=6
	ERML_extlikemin=15
elif Version==4:
	ERML_CUTRAD=14
	ERML_likemin=6
	ERML_extlikemin=15
elif Version==5:
	ERML_CUTRAD=12
	ERML_likemin=6
	ERML_extlikemin=15
elif Version==6 or Version==9:
	ERML_CUTRAD=16
	ERML_likemin=5
	ERML_extlikemin=12
elif Version==7 or Version==10:
	ERML_CUTRAD=14
	ERML_likemin=5
	ERML_extlikemin=12
elif Version==8 or Version==11:
	ERML_CUTRAD=12
	ERML_likemin=5
	ERML_extlikemin=12

if Version in (2,6,7,8):
	ERML_extmin=5
	ERML_extmax=40
	ERML_extlikemin=12
elif Version in (9,10,11):
	ERML_extmin=2
	ERML_extmax=12
	ERML_extlikemin=10
	ERML_MULTRAD=20

if Version==12:
	ERML_CUTRAD=8
	ERML_MULTRAD=20
	ERML_likemin=5
	ERML_extmin=2
	ERML_extmax=12
	ERML_extlikemin=15
elif Version==13: #cluster detection
	ERML_CUTRAD=15
	ERML_MULTRAD=15
	ERML_likemin=5
	ERML_extlikemin=8
	ERML_extmin=2
	ERML_extmax=15
	VerBand = 'C'
	band2use = Bandgroups['C']

if Version==2: Nmlin=2
else: Nmlin=3

inputdir = args.inputdir
if inputdir[-1]=='/': inputdir=inputdir[:-1]
outdir = args.outdir
if outdir=='': outdir = f"{inputdir.replace('_images','_det')}_V{Version:d}{VerBand:s}{'s' if Sixte else ''}"
if not os.path.isdir(outdir): os.mkdir(outdir)
logfile = os.path.join(outdir, f"eRO_V{Version:d}{VerBand:s}{'s' if Sixte else ''}.log")

healpix_id = outdir.split('/')[-1]
outprefix = healpix_id + "_" # ""
print("outprefix",outprefix)
#outprefix = ""

emin_keV = [Bands[i][0] for i in band2use]
emax_keV = [Bands[i][1] for i in band2use]
bandname = [str(i) for i in band2use]
ecf = [Bands[i][2] for i in band2use]

eminstr=' '.join(f"{x*1000:g}" for x in emin_keV)
emaxstr=' '.join(f"{x*1000:g}" for x in emax_keV)
ecfstr= ' '.join(f"{x:g}" for x in ecf)
EvtImgFiles = [os.path.join(inputdir,  f"{outprefix}02{bname}_EvtImg.fits") for bname in bandname]
ExpMapFiles = [os.path.join(inputdir,  f"{outprefix}02{bname}_ExpMap.fits") for bname in bandname]
SrcMapFiles = [os.path.join(outdir,    f"{outprefix}02{bname}_SrcMap.fits") for bname in bandname]
PSFMapFiles = [os.path.join(outdir,    f"{outprefix}02{bname}_PSFMap.fits") for bname in bandname]
BkgMapFiles = {n:[os.path.join(outdir, f"{outprefix}02{bname}_Bg{n:d}Map.fits") for bname in bandname] for n in (1,2,3)}
CheMskFiles = [os.path.join(outdir, f"{outprefix}02{bname}_CheMsk.fits") for bname in bandname]
ApeSenFiles = [os.path.join(outdir, f"{outprefix}02{bname}_ApeSen.fits") for bname in bandname]
UnVigExpMap = [os.path.join(outdir, f"{outprefix}02{bname}_UnvExpMap.fits") for bname in bandname]
DetMask     = [os.path.join(outdir, f"{outprefix}02{bname}_DetMsk.fits") for bname in bandname]

bname = bandname[0]
BoxCats = {n:os.path.join(outdir, f"{outprefix}02{bname}_Bo{n:d}Cat.fits") for n in (1,2,3,4,5)}
MLCats = {n:os.path.join(outdir, f"{outprefix}02{bname}_ML{n}Cat.fits") for n in (1,2,'A')}
MLinCat = {1:BoxCats[Nmlin], 2:os.path.join(outdir, f"{outprefix}02{bname}_ML1Cat.fits"), 'A':os.path.join(outdir, f"{outprefix}02{bname}_ML2Cat.fits")}
ApeCat = os.path.join(outdir, f"{outprefix}02{bname}_ApeCat.fits")
SrcCats = {n:os.path.join(outdir, f"{outprefix}02{bname}_Sc{n:d}Cat.fits") for n in (1,2)}
SrcCat = os.path.join(outdir, f"{outprefix}02{bname}_SrcCat.fits")
CatApe = os.path.join(outdir, f"{outprefix}02{bname}_CatApe.fits") #useless
SkyCovSoft = os.path.join(outdir, f"{outprefix}SkyCov_Soft_{MinDetLike:g}.fits")
SenMapSoft = os.path.join(outdir, f"{outprefix}SenMap_Soft_{MinDetLike:g}.fits")
SkyCovHard = os.path.join(outdir, f"{outprefix}SkyCov_Hard_{MinDetLike:g}.fits")
SenMapHard = os.path.join(outdir, f"{outprefix}SenMap_Hard_{MinDetLike:g}.fits")

cmd_erboxlocal = ["erbox",
	f"images={' '.join(EvtImgFiles)}",
	f"boxlist={BoxCats[1]}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"detmasks={DetMask[0]}",
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

cmd_backmap = dict()
for Nrun in (1,2,3):
	cmd_backmap[Nrun] = []
	for i in range(len(bandname)):
		cmd_backmap[Nrun].append(f"""erbackmap boxlist={BoxCats[Nrun]} bkgimage={BkgMapFiles[Nrun][i]} \
image={EvtImgFiles[i]} expimage={ExpMapFiles[i]} detmask={DetMask[0]} cheesemask={CheMskFiles[i]} \
idband={i+1:d} scut=0.0001 mlmin=6 maxcut=0.5 \
fitmethod=smooth nsplinenodes=36 degree=2 smoothflag=yes \
smoothval=15.  snr=40.0 excesssigma=1000.  nfitrun=3 cheesemaskflag=Y
""")

cmd_erbox = dict()
for Nrun in (2,3):
	cmd_erbox[Nrun] = f"""erbox boxlist="{BoxCats[Nrun]}" \
images="{' '.join(EvtImgFiles)}" expimages="{' '.join(ExpMapFiles)}" bkgimages="{' '.join(BkgMapFiles[Nrun-1])}" \
emin="{eminstr}" emax="{emaxstr}" ecf="{ecfstr}" \
likemin=4 boxsize=4 nruns=2 detmasks="{DetMask[0]}" \
compress_flag=N bkgima_flag=Y expima_flag=Y detmask_flag=Y
"""

cmd_ermldet = dict()
for Nrun in (1,2):
	cmd_ermldet[Nrun] = f"""ermldet boxlist={MLinCat[Nrun]} mllist={MLCats[Nrun]} \
images="{' '.join(EvtImgFiles)}" expimages="{' '.join(ExpMapFiles)}" bkgimages="{' '.join(BkgMapFiles[Nmlin])}" \
emin="{eminstr}" emax="{emaxstr}" ecf="{ecfstr}" \
detmasks="{DetMask[0]}" srcimages="{' '.join(SrcMapFiles)}" \
likemin={ERML_likemin} extlikemin={ERML_extlikemin} cutrad={ERML_CUTRAD} multrad=${ERML_MULTRAD} \
compress_flag=N extmin={ERML_extmin} extmax={ERML_extmax} bkgima_flag=Y expima_flag=Y \
detmask_flag=Y extentmodel=beta thres_flag=N thres_col=like thres_val=30. \
nmaxfit=3 nmulsou=2 fitext_flag=yes srcima_flag=yes \
shapelet_flag=yes photon_flag=yes sixte_flag={'yes' if Sixte else 'no'} extlike_slope="0.0 0.0" \
&& mv {MLCats[Nrun]} {MLCats[Nrun]}.save -v \
&& ftcopy "{MLCats[Nrun]}.save[col *,RADEC_ERR=(RADEC_ERR<=0.0?(sqrt(2.)*max(X_IMA_ERR,Y_IMA_ERR)*4):RADEC_ERR)]" {MLCats[Nrun]} clobber=yes
"""

cmd_catprep = {n:f"catprep infile={MLCats[n]} outfile={SrcCats[n]}" for n in (1,2)}

cmd_apetool=["apetool",
	f"mllist={MLCats['A']}",
	f"apelist={MLCats['A']}",
	f"apelistout={ApeCat}",
	f"images={' '.join(EvtImgFiles)}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"bkgimages={' '.join(BkgMapFiles[Nmlin])}",
	f"psfmaps={' '.join(PSFMapFiles)}",
	f"detmasks={DetMask[0]}",
	f"apesenseimages={' '.join(ApeSenFiles)}",
	f"srcimages={' '.join(SrcMapFiles)}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"eindex=",
	F"eefextract=0.75",
	f"pthresh=4e-6",
	F"cutrad={ERML_CUTRAD}",
	f"psfmapsampling=11",
	F"apexflag=yes",
	f"psfmapflag=yes",
	f"apesenseflag=yes",
	f"shapepsf=yes"]

cmd_apetool16=["apetool",
	f"mllist={MLCats['A'].replace('.fits','_16.fits')}",
	f"apelist={MLCats['A'].replace('.fits','_16.fits')}",
	f"apelistout={ApeCat.replace('.fits','_16.fits')}",
	f"images={' '.join(EvtImgFiles)}",
	f"expimages={' '.join(ExpMapFiles)}",
	f"bkgimages={' '.join(BkgMapFiles[Nmlin])}",
	f"psfmaps={' '.join(PSFMapFiles)}",
	f"detmasks={DetMask[0]}",
	f"apesenseimages={' '.join(ApeSenFiles)}",
	f"srcimages={' '.join(SrcMapFiles)}",
	f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV),
	f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV),
	f"eindex=",
	F"eefextract=0.75",
	f"pthresh=4e-6",
	F"cutrad=16",
	f"psfmapsampling=11",
	F"apexflag=yes",
	f"psfmapflag=yes",
	f"apesenseflag=yes",
	f"shapepsf=yes"]

#cmd_sensmap_s = ["ersensmap", f"sensimage={SenMapSoft}", f"area_table={SkyCovSoft}", f"expimages={' '.join(ExpMapFiles)}", f"detmasks={DetMask[0]}", f"bkgimages={' '.join(BkgMapFiles)}", f"srcimages={' '.join(SrcMapFiles)}", f"emin=" + " ".join(f"{x*1000:g}" for x in emin_keV), f"emax=" + " ".join(f"{x*1000:g}" for x in emax_keV), f"ecf=" + " ".join(f"{x:g}" for x in ecf_soft), f"method=FIT", f"likemin={MinDetLike:g}", f"extflag=no", f"extentmodel=beta", f"extlikemin=15", f"detmask_flag=yes", f"shapelet_flag=no", f"photon_flag=no", f"area_flag=yes" ]

################################################################################
localvars = locals()
#if os.path.isfile(logfile): sys.exit(f'ERROR: {logfile} already exists, please delete it')
Parsdict = dict((k, localvars[k]) for k in ['Version', 'inputdir', 'outdir', 'outprefix', 'cmd_erboxlocal', 'cmd_backmap', 'cmd_erbox', 'cmd_ermldet', 'cmd_catprep' ])
Parsdict['History'] = __doc__.split('\n')
with open(logfile, 'w') as flog:
	json.dump(Parsdict, flog, indent=4)

preclear(BoxCats[1])
subprocess.run(cmd_erboxlocal, check=True)

preclear(BkgMapFiles[1], CheMskFiles)
for cmd in cmd_backmap[1]:
	runcmd(cmd,logfile)

preclear(BoxCats[2])
runcmd(cmd_erbox[2], logfile)

preclear(BkgMapFiles[2], CheMskFiles)
for cmd in cmd_backmap[2]:
	runcmd(cmd, logfile)

preclear(BoxCats[3])
runcmd(cmd_erbox[3],logfile)

preclear(BkgMapFiles[3], CheMskFiles)
for cmd in cmd_backmap[3]: runcmd(cmd,logfile)

assert all(os.path.isfile(bm) for bm in BkgMapFiles[Nmlin])
preclear(SrcMapFiles, MLCats[1])
runcmd(cmd_ermldet[1],logfile)
if not os.path.isfile(MLCats[1]): sys.exit('ERROR: failed')

preclear(SrcCat[1])
runcmd(cmd_catprep[1],logfile)

os.system("chgrp erosim "+outdir+"/*")

sys.exit()

preclear(SrcMapFiles, MLCats[2])
runcmd(cmd_ermldet[2],logfile)
if not os.path.isfile(MLCats[2]): sys.exit('ERROR: failed')

preclear(SrcCat[2])
runcmd(cmd_catprep[2],logfile)

preclear(ApeSenFiles,PSFMapFiles,ApeCat)
subprocess.call(f"cp {MLCats[1]} {MLCats['A']}", shell=True)
subprocess.run(cmd_apetool, check=True)

preclear(ApeSenFiles,PSFMapFiles,ApeCat.replace('.fits','_16.fits'))
subprocess.call(f"cp {MLCats[1]} {MLCats['A'].replace('.fits','_16.fits')}", shell=True)
subprocess.run(cmd_apetool16, check=True)
