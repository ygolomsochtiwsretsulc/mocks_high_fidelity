#!/usr/bin/env python

import os
import astropy.io.fits as fits
import sys
import numpy as np
import glob
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import json
from scipy.stats import norm
from argparse import ArgumentParser



# reading arguments
parser = ArgumentParser()
parser.add_argument('-o','--outroot',help='output root file (default=None)',type=str,default=None,metavar='OUTROOT')
parser.add_argument(     '--dz',help='zspec binning for GMM (default=0.05)',type=float,default=0.05,metavar='DZ')
parser.add_argument(     '--zmin',help='zmin (default=0.0)',type=float,default=0.0,metavar='ZMIN')
parser.add_argument(     '--zmax',help='zmax (default=2.0)',type=float,default=2.0,metavar='ZMAX')
parser.add_argument('-n','--ngauss',help='number of gaussians (default=5)',type=int,default=5,metavar='NGAUSS')
args        = parser.parse_args()

# selections
sels = ['elg']

# zspec grid
zgrid = np.array([args.zmin+i*args.dz for i in range(10000) if args.zmin+i*args.dz<=args.zmax])
zcents= 0.5*(zgrid[1:]+zgrid[:-1])
nbins = len(zgrid)-1

'''
# kids files
fns = glob.glob('/hpcstorage/raichoor/kids_dr4/KiDS_DR4.*_ugriZYJHKs_cat.fits')
#fns = fns[:100]

# iselg sel
def get_issel(name,g,gr,ri):
	if (name=='elg'):
		keep  = (g>21.0) & (g<23.2)
		keep &= (0.5-2.5*gr<ri) & (ri<2.5-2.5*gr) & (0.4*gr+0.3<ri) & (ri<0.4*gr+0.9)
	return keep


keys  = ['id','kids_tile','alpha_j2000','delta_j2000',
			'rtot','rfib','gr','ri','iz','zy','yj','jh','hks','zphot']
# initialising
mydict   = {}
for key in keys + sels:
	mydict[key] = []


for fn in fns:
	print 'reading ', fn
	hdu  = fits.open(fn)
	data = hdu[1].data
	## clean
	keep = (data['mask']==0)
	data = data[keep]
	## rfib: pixel scale=0.2 arcsec/pix (TBC!); so 4 pix radius = 0.8 arcsec radius; 4MOST fibres are 0.7 arcsec radius
	## gtot = gr + rtot
	## colour_gaap are dereddened
	tmpdict = {}
	tmpdict['rtot'] = data['mag_auto']  -data['extinction_r'] 
	tmpdict['rfib'] = data['mag_aper_4']-data['extinction_r']
	tmpdict['zphot']= data['z_b']
	tmpdict['zphot'][data['sg_flag']==0] = 0
	#
	keep  = np.zeros(len(data),dtype=bool)
	for sel in sels:
		tmpdict[sel] = get_issel(sel,
								data['colour_gaap_g_r']+tmpdict['rtot'],
								data['colour_gaap_g_r'],
								data['colour_gaap_r_i'])
		keep |= tmpdict[sel]
	## cutting
	data = data[keep]
	for key in tmpdict.keys():
		tmpdict[key] = tmpdict[key][keep]
	## appending
	for key in keys+sels:
		if (key in ['rtot','rfib','zphot'] + sels):
			mydict[key] += tmpdict[key].tolist()
		elif (key in ['gr','ri','iz','zy','yj','jh','hks']):
			mydict[key] += data['colour_gaap_'+key[0]+'_'+key[1:]].tolist()
		else:
			mydict[key] += data[key].tolist()

# writing fits
collist = []
for key in keys:
	# if the columns directly comes from kids, we take the format from the last read hdu
	if (key.upper() in hdu[1].columns.names):
		fmt = hdu[1].columns[key.upper()].format
	else: # float
		fmt = 'E'
	collist.append(fits.Column(name=key,format=fmt,array=mydict[key]))
for sel in sels:
	collist.append(fits.Column(name='is'+sel,format='L',array=mydict[sel]))
hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
hdu.writeto(args.outroot+'.fits',overwrite=True)
'''


# gaussian mixture
# http://www.astroml.org/book_figures/chapter4/fig_GMM_1D.html
# https://scikit-learn.org/stable/auto_examples/mixture/plot_gmm_selection.html#sphx-glr-auto-examples-mixture-plot-gmm-selection-py
hdu = fits.open(args.outroot+'.fits')
data= hdu[1].data
##
mydict = {}
mydict['dz'],mydict['zmin'],mydict['zmax'] = args.dz,args.zmin,args.zmax
fn     = open(args.outroot+'.nzgmm.json','w')
for sel in sels:
	print sel
	#
	bic,lowest_bic     = [],np.infty
	n_components_range = range(1, args.ngauss)
	#
	z = data['zphot'][data['is'+sel]]
	print sel, len(z), 'obj.'
	# gmm
	X = z.reshape((len(z),1))
	for n_components in range(1,args.ngauss):
		gmm = GaussianMixture(n_components=n_components,covariance_type='diag')
		gmm.fit(X)
		bic.append(gmm.bic(X))
		if bic[-1] < lowest_bic:
			lowest_bic = bic[-1]
			best_gmm = gmm
	bic = np.array(bic)
	clf = best_gmm
	## rounding
	ps   = np.round(clf.weights_,3).tolist()
	mus  = np.round(clf.means_.flatten(),3).tolist()
	sds  = np.round(np.sqrt(clf.covariances_.flatten()),3).tolist()
	mydict[sel] = {}
	mydict[sel]['p'],mydict[sel]['mu'],mydict[sel]['sd'] = ps,mus,sds
	## plotting
	fig,ax = plt.subplots()
	nzraw = np.array([((z>=zgrid[i]) & (z<zgrid[i+1])).sum() for i in range(nbins)])
	x,y    = zcents,nzraw
	ax.bar(x,y,y*0+args.dz,color='k',edgecolor='k',alpha=0.3,label='data')
	dx   = 0.001
	x    = np.arange(args.zmin,args.zmax+dx,dx)
	yfit_indiv = np.array([p*norm.pdf(x,mu,sd) for p,mu,sd in zip(ps,mus,sds)])
	tmpx,tmpy  = 0.55,0.8
	ax.text(tmpx,tmpy+0.05,'dz='+str(args.dz),transform=ax.transAxes)
	for p,mu,sd in zip(ps,mus,sds):
	    ax.text(tmpx,tmpy,'p='+str(p)+'  mu='+str(mu)+'  sd='+str(sd),transform=ax.transAxes)
	    tmpy -= 0.05
	yfit = np.sum(yfit_indiv, axis = 0)
	# np.trapz(yfit,x=x) = 1
	yfit *= np.sum(y)*args.dz
	ax.plot(x, yfit,color='r',label='GMM fit')
	ax.grid(True)
	ax.set_xlabel('kids zphot')
	ax.set_xlim(args.zmin,args.zmax)
	ax.set_title(sel+' ('+str(len(z))+' obj.)')
	ax.legend()
	plt.savefig(args.outroot+'.nzgmm.'+sel+'.png',bbox_inches='tight')
	plt.close()

#
json.dump(mydict,fn,cls=json.JSONEncoder)
fn.close()

