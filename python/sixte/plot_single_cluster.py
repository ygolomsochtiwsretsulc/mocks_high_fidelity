"""
cp -r /data40s/erosim/eRASS/eRASS8_cluster_MD40/figures $GIT_AGN_MOCK/figures/eRASS8_cluster_MD40/
"""
import os, sys, glob
from scipy.interpolate import interp1d
import astropy.io.fits as fits
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import norm
from astropy.table import Table, Column
from astropy_healpix import healpy
import numpy as n

env = 'MD40'
fig_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'figures', env, 'clusters')
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)

# simulation setup
if env[:2] == "MD" : # env == "MD04" or env == "MD40" or env == "MD10" or env == "MD25"
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env[:4] == "UNIT" : # == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
    from astropy.cosmology import FlatLambdaCDM
    import astropy.units as u
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT


import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

cluster_file = os.path.join(os.environ[env], env+'_eRO_CLU.fit')
clu = Table.read(cluster_file)

event_file = "/data40s/erosim/eRASS/eRASS8_cluster_MD40/simulated_photons_ccdA.fits"
evt_all = fits.open(event_file)[1].data

simput_dir = os.path.join(os.environ[env], 'cat_CLU_SIMPUT')

HEALPIX_8 = healpy.ang2pix(8, n.pi/2. - clu['DEC']*n.pi/180. , clu['RA']*n.pi/180. , nest=True)
hpx_val = 000
in_simput = clu[(HEALPIX_8 == hpx_val)]
simput_file = os.path.join(simput_dir, str(hpx_val).zfill(6)+'.fit')
spt = Table.read(simput_file)

event_files = n.array(glob.glob(os.path.join( '/data40s/erosim/eRASS/eRASS8_cluster_MD40', str(hpx_val).zfill(3), 'erass_ccd?_evt.fits' )))

def get_events(src_id):
	#src_id = simput_entry['SRC_ID']
	all_ev, ras, decs = [], [], []
	for evf in event_files:
		hdu = fits.open(evf)
		s1 = (hdu[1].data['SRC_ID'].T[0]==src_id)
		all_ev.append(hdu[1].data['SIGNAL'][s1])
		ras.append(hdu[1].data['RA'][s1])
		decs.append(hdu[1].data['DEC'][s1])

	all_ev = n.hstack((all_ev))
	return all_ev, n.hstack((ras)), n.hstack((decs))

# select cluster id : 
#id_clu = 456
middle = int(len(spt['FLUX'])*2)
middle2 = int(len(spt['FLUX'])*3/4)
ids_loop = n.hstack(( n.where(in_simput['HALO_pid']>0)[0], n.argsort(spt['FLUX'])[:5] , n.argsort(spt['FLUX'])[middle:middle+10], n.argsort(spt['FLUX'])[middle2:middle2+10], n.argsort(spt['FLUX'])[-10:] ))
for id_clu in ids_loop:
	print(id_clu)
	catalog_entry = in_simput[id_clu]
	simput_entry = spt[id_clu]

	path_2_image = os.path.join( simput_dir, simput_entry['IMAGE'][:-7] )
	path_2_spec = os.path.join( simput_dir, simput_entry['SPECTRUM'][:-19] )
	image = fits.open(path_2_image)[0].data
	spectrum = fits.open(path_2_spec)[1].data

	evt = evt_all[(evt_all['SRC_ID_1'] == simput_entry['SRC_ID'])]
	all_evt_kev, ras_evt, decs_evt = get_events(simput_entry['SRC_ID'])
	basename = os.path.join(fig_dir, simput_entry['IMAGE'][:-12].split('/')[-1])

	# panel 1 
	# plot image 

	p.figure(1, (5,5))
	p.imshow(image)
	#p.imshow(n.log10(image))
	p.colorbar()
	p.title('surface brightness')
	p.savefig(basename + '-image.png')
	p.clf()

	# panel 2 
	# plot events
	# energy histogram along with spectrum 

	bins = 10**n.arange(-1, 2., 0.05)
	p.figure(1, (5,5))
	p.axes([0.18, 0.16, 0.75, 0.75])
	p.plot( spectrum['ENERGY'][0], spectrum['FLUXDENSITY'][0]/n.median(spectrum['FLUXDENSITY'][0]), label=os.path.basename(path_2_spec)[17:] )
	p.hist(all_evt_kev, bins=bins, normed=True, histtype='step')
	p.xscale('log')
	p.xlim((bins[0], bins[-1]))
	#p.ylim((catalog_entry['DEC']-5/60., catalog_entry['DEC']+5/60.))
	p.xlabel('E [keV]')
	p.ylabel('Normed counts [d log10(keV)=0.05] ')
	p.legend(loc=0)
	p.savefig(basename + '-spectrum.png')
	p.clf()

	# ra-dec distribution 
	# 1x, 2x r500c and rvir into circles / ellipses
	rvir_deg = (catalog_entry['HALO_Rvir'] *  cosmo.arcsec_per_kpc_proper (catalog_entry['redshift_R']) / 3600.).value
	theta = n.linspace(0, 2*n.pi, 100)

	p.figure(1, (5,5))
	p.axes([0.18, 0.16, 0.75, 0.75])
	p.xlabel('R.A. [deg]')
	p.ylabel('Dec. [deg]')
	x1 = rvir_deg * n.cos(theta) + catalog_entry['RA']
	x2 = rvir_deg * n.sin(theta) + catalog_entry['DEC']
	p.plot(x1,x2,'k--', label='1 r vir')

	x1 = 2 * rvir_deg * n.cos(theta) + catalog_entry['RA']
	x2 = 2 * rvir_deg * n.sin(theta) + catalog_entry['DEC']
	p.plot(x1,x2,'b--', label='2 r vir')

	p.plot(in_simput['RA'], in_simput['DEC'], 'k,')
	p.plot(ras_evt, decs_evt, 'b+' , label=str(len(ras_evt))+' events', alpha=0.5)
	p.plot(catalog_entry['RA'], catalog_entry['DEC'], 'ro', label='' )
	#p.plot(ras_evt, decs_evt, 'bs',color='black', markerfacecolor='none', markeredgecolor='black', markeredgewidth=1)
	p.xlim((catalog_entry['RA']-3.*rvir_deg, catalog_entry['RA']+3.*rvir_deg))
	p.ylim((catalog_entry['DEC']-3.*rvir_deg, catalog_entry['DEC']+3.*rvir_deg))
	p.title('M500c='+str( n.round(n.log10(catalog_entry['HALO_M500c']),2) ) + ', z=' + str( n.round(catalog_entry['redshift_R'],2) ))
	p.legend(loc=0)
	p.savefig(basename + '-ra-dec.png')
	p.clf()

# contour plot with number density of counts ? to reflec the image profile ?

os.system("cp /home/comparat/software/linux/lss_mock_dev/figures/MD40/clusters/*.png /home/comparat/wwwDir/stuff/ ")

# panel 3
# plot density profile



