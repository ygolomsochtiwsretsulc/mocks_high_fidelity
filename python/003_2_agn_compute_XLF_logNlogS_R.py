"""
What it does
------------

Computes and tabulates the logNlogS, logNlogR, 

Creates figures with the tabulate data

References
----------

Command to run
--------------

python3 003_2_agn_compute_XLF_logNlogS_R.py environmentVAR ftyp

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

ftyp: 'all' or 'sat', type of AGN populating a central halo or satellite haloes.

Dependencies
------------

import time, os, sys, numpy, scipy, astropy, h5py, matplotlib

"""
import time
t0 = time.time()

import glob, os, sys

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p

import astropy.io.fits as fits
from astropy_healpix import healpy
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
from scipy.interpolate import interp1d

#import h5py
import numpy as n

print('TABULATES XLF and logNlogS')
print('------------------------------------------------')
print('------------------------------------------------')

#import astropy.io.fits as fits
# import all pathes

env = sys.argv[1] # 'MD04'
print(env)


# simulation setup
if env[:2] == "MD" :
    cosmoMD = FlatLambdaCDM(
        H0=67.77 * u.km / u.s / u.Mpc,
        Om0=0.307115)  # , Ob0=0.048206)
    h = 0.6777
    L_box = 1000.0 / h
    cosmo = cosmoMD
if env[:4] == "UNIT_fA1_DIR" :
    cosmoUNIT = FlatLambdaCDM(H0=67.74 * u.km / u.s / u.Mpc, Om0=0.308900)
    h = 0.6774
    L_box = 1000.0 / h
    cosmo = cosmoUNIT

#if env == "MD10" or env == "UNIT_fA1_DIR" or env == "UNIT_fA1i_DIR" or env == "UNIT_fA2_DIR":
z_maximum = 6.0
if env == "MD04":
    z_maximum = 0.45


agn_catalog_list = sorted( n.array( glob.glob( os.path.join( os.environ[env], 'cat_AGN_all', '*.fit' ) ) ) )
agn_catalog_list.sort()

galaxy_catalog_list = n.array( glob.glob( os.path.join( os.environ[env], 'cat_GALAXY_all', '*.fit' ) ) )
galaxy_catalog_list.sort()

agn_data_dir = os.path.join( os.environ['GIT_AGN_MOCK'], 'data', 'agn' )
fig_dir      = os.path.join( os.environ['GIT_AGN_MOCK'], 'figures', env, 'agn' )
logNlogS_dir = os.path.join( fig_dir, 'logNlogS' )
logNlogR_dir = os.path.join( fig_dir, 'logNlogR' )
XLF_dir      = os.path.join( fig_dir, 'XLF')
DC_dir       = os.path.join( fig_dir, 'duty_cycle')
LSAR_dir     = os.path.join( fig_dir, 'LSAR')
if os.path.isdir(fig_dir) == False:
    os.system('mkdir -p ' + fig_dir)
if os.path.isdir(logNlogS_dir) == False:
    os.system('mkdir -p ' + logNlogS_dir)
if os.path.isdir(logNlogR_dir) == False:
    os.system('mkdir -p ' + logNlogR_dir)
if os.path.isdir(XLF_dir) == False:
    os.system('mkdir -p ' + XLF_dir)
if os.path.isdir(DC_dir) == False:
    os.system('mkdir -p ' + DC_dir)
if os.path.isdir(LSAR_dir) == False:
    os.system('mkdir -p ' + LSAR_dir)

area = healpy.nside2pixarea(8, degrees=True)

f_duty = interp1d(n.array([0., 0.75, 2., 3.5, 10.1]), n.array([0.1, 0.2, 0.3, 0.3, 0.3]))

# redshift binning 
DZ = 0.2
all_zs = n.arange(0., z_maximum, DZ)

# logN logS: X-ray flux binning
fx_bins = n.arange(-20, -8., 0.2)
x_fx = fx_bins[:-1] + 0.1

# logN logR: r magnitude binning
dlog_magr = 0.1
mag_r_bins = n.arange(8, 40 + 2 * dlog_magr, dlog_magr)
x_mag_r = mag_r_bins[:-1]+dlog_magr/2.

# LF binning
dlog_lx = 0.05
lx_bins = n.arange(36, 48, dlog_lx)
x_lx = lx_bins[:-1] + dlog_lx / 2.

# LSAR histogram
dlog_LSAR = 0.1
LSAR_bins = n.arange(30, 38, dlog_LSAR)
x_LSAR = LSAR_bins[:-1] + dlog_LSAR / 2.

# SMF and duty cycle
dlogM = 0.1
bins_SMF = n.arange(8, 13, dlogM)
x_SMF = (bins_SMF[1:] + bins_SMF[:-1]) * 0.5

# galaxy stellar mass function of Ilbert 2013
# mass bins 
mbins = n.arange(8,12.5,0.25)
path_ilbert13_SMF = os.path.join( os.environ['GIT_AGN_MOCK'], 'data','LF_SMF', "ilbert_2013_mass_function_params.txt")
smf_ilbert13 = lambda M, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s : ( phi_1s * (M/M_star) ** alpha_1s + phi_2s * (M/M_star) ** alpha_2s ) * n.e ** (-M/M_star) * (M/ M_star)
smf_ilbert_zmin, smf_ilbert_zmax, N, M_comp, M_star, phi_1s, alpha_1s, phi_2s, alpha_2s, log_rho_s = n.loadtxt(path_ilbert13_SMF, unpack=True)
smf_ilbert_fun = n.array([
lambda mass : smf_ilbert13( mass , 10**M_star[0], phi_1s[0]*10**(-3), alpha_1s[0], phi_2s[0]*10**(-3), alpha_2s[0] )
, lambda mass : smf_ilbert13( mass , 10**M_star[1], phi_1s[1]*10**(-3), alpha_1s[1], phi_2s[1]*10**(-3), alpha_2s[1] )
, lambda mass : smf_ilbert13( mass , 10**M_star[2], phi_1s[2]*10**(-3), alpha_1s[2], phi_2s[2]*10**(-3), alpha_2s[2] )
, lambda mass : smf_ilbert13( mass , 10**M_star[3], phi_1s[3]*10**(-3), alpha_1s[3], phi_2s[3]*10**(-3), alpha_2s[3] )
, lambda mass : smf_ilbert13( mass , 10**M_star[4], phi_1s[4]*10**(-3), alpha_1s[4], phi_2s[4]*10**(-3), alpha_2s[4] )
, lambda mass : smf_ilbert13( mass , 10**M_star[5], phi_1s[5]*10**(-3), alpha_1s[5], phi_2s[5]*10**(-3), alpha_2s[5] )
, lambda mass : smf_ilbert13( mass , 10**M_star[6], phi_1s[6]*10**(-3), alpha_1s[6], phi_2s[6]*10**(-3), alpha_2s[6] )
, lambda mass : smf_ilbert13( mass , 10**M_star[7], phi_1s[7]*10**(-3), alpha_1s[7], phi_2s[7]*10**(-3), alpha_2s[7] )
])

smf_ilbert_name = n.array([ "Il13 "+str(zmin)+"<z<"+str(zmax) for zmin, zmax in zip(smf_ilbert_zmin,smf_ilbert_zmax) ])

def get_lognlogs(hpx_id='000061'):
    ff = os.path.join(logNlogS_dir, 'logNlogS_soft_' + hpx_id + '.ascii')
    # print(ff)
    outout = n.loadtxt(ff, unpack=True)
    #print(outout, outout.shape)
    # data=n.array(data)
    #N_pixels = len(data)
    NN = outout / area
    #itp = interp1d(x_fx, n.log10(NN))
    return x_fx, n.log10(NN) # itp

for path_2_eRO_catalog, path_2_GAL_catalog in zip(agn_catalog_list[::100], galaxy_catalog_list[:100]):
	print(path_2_eRO_catalog, path_2_GAL_catalog, time.time()-t0)
	str_healpix_id = os.path.basename(path_2_eRO_catalog)[:-4]
	hd = fits.open(path_2_eRO_catalog)
	hg = fits.open(path_2_GAL_catalog)
	# AGN data
	z = hd[1].data['redshift_R']
	lx = hd[1].data['LX_hard']
	logNH = hd[1].data['logNH']
	fx = hd[1].data['FX_soft']
	lx_0520 = hd[1].data['LX_soft']
	logm = hd[1].data['galaxy_SMHMR_mass']
	lsar = lx - logm
	mag_r = hd[1].data['SDSS_r_AB']
	g_lat = hd[1].data['g_lat']
	# galaxy data
	logm_gal = hg[1].data['galaxy_SMHMR_mass']
	z_gal = hg[1].data['redshift_R']
	#n_agn = len(z)

	print(" TABULATE LOG N LOG S, fx ", time.time()-t0)
	def get_lognlogs_replicas(fx, area):
		log_f_05_20 = n.log10(fx[fx > 0])
		out = n.histogram(log_f_05_20, bins = fx_bins )
		# cumulative number density per square degrees
		N_out = n.cumsum(out[0][::-1])[::-1]
		#c_out = N_out / area
		#c_out_up = (1 + N_out**(-0.5)) * c_out
		#c_out_low = (1 - N_out**(-0.5)) * c_out
		#c_err = (n.log10(c_out_up) - n.log10(c_out_low)) / 2.
		return N_out  # , c_err

	#Xgal = (abs(g_lat)>20)
	N_out = get_lognlogs_replicas(fx, area=area)
	n.savetxt( os.path.join( logNlogS_dir, 'logNlogS_soft_' + str_healpix_id + '.ascii'), N_out ) 

	p.figure(1, (6, 6))
	p.plot(x_fx, n.log10(N_out/area), rasterized=True, lw=2, ls='solid', label='MD10 agn') 
	# Georgakakis 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_VS"],
		'data',
		'logNlogS',
		'logNlogS_Georgakakis_08_AGN.data')
	x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
	p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='g', label='G08')

	# Merloni 2012
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_VS"],
		'data',
		'logNlogS',
		'logNlogS_Merloni_12_AGN.data')
	x_data, y_data = n.loadtxt(path_2_logNlogS_data, unpack=True)
	p.plot(n.log10(x_data), n.log10(y_data), lw=3, ls='dotted', color='r', label='M12')

	# Mateos 2008
	path_2_logNlogS_data = os.path.join(
		os.environ["GIT_VS"],
		'data/logNlogS/logNlogS_Mateos_08_AGN.data')
	x_data, y_data, err = n.loadtxt(path_2_logNlogS_data, unpack=True)
	p.plot(x_data, n.log10(y_data), lw=3, ls='dotted', color='b', label='M08')

	p.xlabel('log(F_X[0.5-2 keV])')
	p.ylabel('log(>F_X) [/deg2]')
	p.legend(frameon=False, loc=0)
	# p.yscale('log')
	p.xlim((-18, -11.5))
	p.ylim((-2, 4.2))
	# p.title('Mocks')
	p.grid()
	p.savefig(os.path.join(logNlogS_dir, "logN_logS_AGN" + str_healpix_id + ".png"))
	p.clf()
	print(" TABULATE LOG N LOG R, mag_r ", time.time()-t0)
	def get_lognlogr(mag_r):
		# cumulative number density per square degrees
		out = n.histogram(mag_r, bins = mag_r_bins)
		N_out = n.cumsum(out[0])
		return N_out  

	N_out = get_lognlogr(mag_r)
	# DATA_R = n.transpose(out)#[(out[1]>0) & (out[2]!=n.inf)]
	n.savetxt( os.path.join( logNlogR_dir, 'logNlogR_optical_' + str_healpix_id + '.ascii' ) , N_out )

	p.figure(1, (6,6))
	#p.axhline(n.log10(300), ls='dashed')
	ff_lims = [7e-17, 8e-16, 6e-15, 2e-15, 2e-14, 1e-14]
	ff_lim_labs = ["7e-17", "8e-16", "6e-15", "2e-15", "2e-14", "1e-14"]
	for ff_lim, ff_lab in zip(ff_lims, ff_lim_labs) : 
		NN = get_lognlogr( mag_r[(fx>ff_lim)] )
		p.plot(x_mag_r, n.log10(NN/area), rasterized = True, label = ff_lab  , lw=2, ls='dashed')

	#is_agn = (log_f_05_20>n.log10(7e-17)) & (A_xmm) & (ok) & ( t1 | t2 | ell )
	#x_out_0, c_out_0, c_err_0, c_out = get_lognlogr(is_agn, area_xmm)
	x_data, y_data = n.loadtxt(os.path.join(  os.environ['GIT_AGN_MOCK'], 'data', 'logNlogR', 'logNlogR_S82X_FX_gt_7e-17.ascii'), unpack = True )
	err_pc = (6*10**y_data)**(-0.5)
	p.fill_between(x_data, y1=y_data*(1-err_pc), y2=y_data*(1+err_pc), color='k', label='S82X', alpha=0.3 )

	#is_agn = (log_f_05_20>n.log10(8e-16)) & (A_chandra)  & (ok) & ( t1 | t2 | ell )
	#x_out_0, c_out_0, c_err_0, c_out = get_lognlogr(is_agn, area_chandra)
	#x_data, y_data = n.loadtxt(os.path.join(plotDir, 'logNlogR_S82X_FX_gt_8e-16.ascii'), unpack = True )
	#p.plot(x_data, y_data, lw=3, ls='dotted', color='r', label='S82X FX>8e-16'  )

	##is_agn = (log_f_05_20>n.log10(2e-15)) & (A_10_13) & (ok) & ( t1 | t2 | ell )
	##x_out_0, c_out_0, c_err_0, c_out = get_lognlogr(is_agn, area_10_13)
	#x_data, y_data = n.loadtxt(os.path.join(plotDir, 'logNlogR_S82X_FX_gt_2e-15.ascii'), unpack = True )
	#p.plot(x_data, y_data, lw=3, ls='dotted', color='m', label='S82X FX>2e-15'  )
	p.axvline(23.5, color='b', ls='dotted')
	p.xlabel('magnitude r (SDSS)')
	p.ylabel('log(<m) [/deg2]')
	p.legend(frameon=False, loc=0)
	#p.yscale('log')
	p.xlim((12, 25))
	p.ylim((-3, 3))
	#p.title('Mocks')
	p.grid()
	p.savefig(os.path.join(logNlogR_dir, "logN_logR_AGN" + str_healpix_id + ".png"))
	p.clf()


	for id_z in n.arange(len(all_zs)):
		zmin = all_zs[id_z]
		zmax = all_zs[id_z] + DZ
		z_mean = 0.5 * (zmin + zmax)
		z_selection = (z >= zmin) & (z < zmax)
		z_selection_GAL = (z_gal >= zmin) & (z_gal < zmax)
		baseName = str_healpix_id + '_' + \
			str(n.round(zmin, 2)) + '_z_' + str(n.round(zmax, 2))

		print(zmin, '<z<', zmax, ', t=', time.time()-t0)
		vol = (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) * n.pi * area / 129600.
		DL_mean_z = (cosmo.luminosity_distance(z_mean).to(u.cm)).value
		#print('volume', vol, 'Mpc3')

		# hard X-ray luminosity function Aird 2015)
		def kz_h(z): return 10**(-4.03 - 0.19 * (1 + z))
		def Ls_h(z): return 10**(44.84 - n.log10(((1 + 2.0) / (1 + z))** 3.87 + ((1 + 2.0) / (1 + z))**(-2.12)))
		def phi_h(L, z): return kz_h(z) / ((L / Ls_h(z))**0.48 + (L / Ls_h(z))**2.27)

		# Selections for the histogram
		NH20 = (z_selection) & (logNH < 22)
		NH22 = (z_selection) & (logNH >= 22) & (logNH < 24)
		NH24 = (z_selection) & (logNH >= 24)

		N_nh20 = n.histogram(lx[NH20], lx_bins)[0] / vol / dlog_lx
		N_nh22 = n.histogram(lx[NH22], lx_bins)[0] / vol / dlog_lx
		N_nh24 = n.histogram(lx[NH24], lx_bins)[0] / vol / dlog_lx

		nhar = n.histogram(lx[z_selection], lx_bins)[0] / vol / dlog_lx
		nharN = n.histogram(lx[z_selection], lx_bins)[0]  # /vol/dlog_lx

		p.figure(1, (6, 6))
		p.axes([0.2, 0.18, 0.75, 0.75])
		#sel = ((zmin+zmax)*0.5>z0)&((zmin+zmax)*0.5<z1)
		# if len(sel.nonzero()[0])>0:
		# if len(list(set(z_center[sel])))>=2:
		#sel = (sel) & (z_center==z_center[sel][0])
		#p.plot(Lx_c[sel], phi[sel], label='Ha05 '+str(n.round(z_center[sel][0],2)), ls='dashed')

		# Aird 2015
		z_mean = (zmin + zmax) * 0.5 * n.ones_like(x_lx)
		# mock
		p.plot(x_lx, phi_h(10**x_lx, z_mean), c='cyan', ls='dashed', lw=2, label='Ai15')  # Aird 2-10 keV LADE')
		p.fill_between(x_lx,
						y1=(nhar) * (1 - (nharN)**(-0.5)),
						y2=(nhar) * (1 + (nharN)**(-0.5)),
						color='g',
						alpha=0.7,
						label='Mock',
						lw=2)  # 2-10  keV')
		p.plot(x_lx, N_nh20, label='nH<22', ls='dotted', lw=2)
		p.plot(x_lx, N_nh22, label='22<nH<24', ls='dotted', lw=2)
		p.plot(x_lx, N_nh24, label='24<nH', ls='dotted', lw=2)
		p.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
		p.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
		p.legend(loc=3, fontsize=14)
		p.yscale('log')
		p.xlim((37., 46.5))
		p.ylim((1 / (2 * vol), 1e-2))
		p.title(str(n.round(zmin, 2)) + "<z<" + str(n.round(zmax, 2)))
		p.grid()
		p.savefig(os.path.join(XLF_dir, "XLF_soft_" + baseName + ".png"))
		p.clf()
		DATA_XLF = n.transpose([x_lx, nhar, (nhar)*(1-(nharN)**(-0.5)), (nhar)*(1+(nharN)**(-0.5)), N_nh20, N_nh22, N_nh24, phi_h(10**x_lx,z_mean)])
		n.savetxt(os.path.join(XLF_dir, 'XLF_soft_'+baseName+'.ascii'), DATA_XLF  )
		# XLF_ratio_
		p.figure(1, (6, 6))
		p.axes([0.18, 0.18, 0.75, 0.75])
		p.plot(x_lx, nhar / phi_h(10**x_lx, z_mean), label='2-10 keV')
		p.fill_between(x_lx,
						y1=1 - (nharN)**(-0.5),
						y2=1 + (nharN)**(-0.5),
						color='green',
						alpha=0.3,
						label='mock ERR')
		p.xlabel(r'$\log_{10}(L^{2-10\, keV}_X/[erg/s])$')
		p.ylabel(r'mock/model')
		p.legend(frameon=False, loc=0)
		p.xlim((37., 46.5))
		p.ylim((0.7, 1.3))
		p.title(str(n.round(zmin, 2)) + "<z<" + str(n.round(zmax, 2)))
		p.grid()
		p.savefig(os.path.join(XLF_dir, "XLF_ratio_" + baseName + ".png"))
		p.clf()
		
		# LSAR histogram
		nall = n.histogram(lsar[z_selection], LSAR_bins)[0] / vol / dlog_LSAR
		nallN = n.histogram(lsar[z_selection], LSAR_bins)[0]

		zsel = (logm >= 12) & (z_selection)
		nall_12 = n.histogram(lsar[zsel], LSAR_bins)[0] / vol / dlog_LSAR
		nallN_12 = n.histogram(lsar[zsel], LSAR_bins)[0]

		zsel = (logm >= 11) & (logm < 12) & (z_selection)
		nall_11 = n.histogram(lsar[zsel], LSAR_bins)[0] / vol / dlog_LSAR
		nallN_11 = n.histogram(lsar[zsel], LSAR_bins)[0]

		zsel = (logm >= 10) & (logm < 11) & (z_selection)
		nall_10 = n.histogram(lsar[zsel], LSAR_bins)[0] / vol / dlog_LSAR
		nallN_10 = n.histogram(lsar[zsel], LSAR_bins)[0]

		zsel = (logm >= 9) & (logm < 10) & (z_selection)
		nall_9 = n.histogram(lsar[zsel], LSAR_bins)[0] / vol / dlog_LSAR
		nallN_9 = n.histogram(lsar[zsel], LSAR_bins)[0]

		p.figure(1, (6, 6))
		p.axes([0.16, 0.15, 0.8, 0.8])

		if z_mean[0]<0.5:
			# Ge17
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z025.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, label='G17', alpha=0.5,color='grey')
			# A18
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_010z050_095M100.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, label='A18', alpha=0.5,color='green')
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_010z050_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, label='A18', alpha=0.5,color='red')
			# We17
			#edd_r = 10**n.arange(-6,0,0.1)
			#p.plot(n.arange(-6,0,0.1)+34, xi_LSAR(edd_r)/dlogf, ls='dashed', lw=2, label='We17 z=0.1' )
			
		if z_mean[0]>=0.5 and z_mean[0]<1. :
			#G17
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z075.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
			#A18
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_050z100_095M100.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,alpha=0.5,color='green')# label='A18, 9.5<M<10', 
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_050z100_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,alpha=0.5,color='red')# label='A18, 10<M<10.5', 
			
		if z_mean[0]>=1. and z_mean[0]<1.5 :
			#G17
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z125.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')# label='G17, M>8',
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_100z150_095M100.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,alpha=0.5,color='green')# label='A18, 9.5<M<10', 
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_100z150_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#label='A18, 10<M<10.5', 
			
		if z_mean[0]>=1.5 and z_mean[0]<2. :
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z175.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_150z200_095M100.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='green')#, label='A18, 9.5<M<10'
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_150z200_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#, label='A18, 10<M<10.5'
			
		if z_mean[0]>=2. and z_mean[0]<2.5 :
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z225.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_200z250_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#, label='A18, 10<M<10.5'
			
		if z_mean[0]>=2.5 and z_mean[0]<3. :
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z275.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_250z300_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm, alpha=0.5,color='red')#, label='A18, 10<M<10.5'
			
		if z_mean[0]>=3 :
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_G17_z350.ascii'), unpack=True)
			fun = interp1d(x, 10**(0.5*(y_min+y_max)) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=10**(y_min)/nrm,  y2=10**(y_max)/nrm, alpha=0.5,color='grey')#, label='G17, M>8'
			x, y_min, y_max = n.loadtxt( os.path.join( agn_data_dir,  'lsar_hist_A17_300z400_100M105.ascii'), unpack=True)
			fun = interp1d(x, 0.5*(y_min+y_max) )
			nrm = quad(fun,x.min(), x.max())[0]
			p.fill_between(x+34, y1=(y_min)/nrm,  y2=(y_max)/nrm,  alpha=0.5,color='red')#label='A18, 10<M<10.5',


		fun = interp1d(x_LSAR, nall)
		nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
		p.plot(x_LSAR, nall / nrm, 'k', lw=3)  # , label='mock all'

		fun = interp1d(x_LSAR, nall_9)
		nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
		p.plot(x_LSAR, nall_9 / nrm, 'g', lw=2)  # , label='9-10'

		fun = interp1d(x_LSAR, nall_10)
		nrm = quad(fun, x_LSAR.min(), x_LSAR.max())[0]
		p.plot(x_LSAR, nall_10 / nrm, 'r', lw=2)  # , label='10-11'

		p.xlabel(r'$\log_{10}(\lambda_{SAR})$')
		p.ylabel(r'probability distribution function')
		#p.legend(frameon=False, loc=3)
		p.yscale('log')
		p.xlim((30., 35.5))
		p.ylim((1e-4, 4))
		p.title('Specific accretion rate, ' + str(n.round(zmin, 2)) + r"<z<" + str(n.round(zmax, 2)))
		p.grid()
		p.savefig(os.path.join(LSAR_dir, "LSAR_hist_" + baseName + ".png"))
		p.clf()
		
		# DUTY CYCLE
		# print("duty_cycle_AGN")
		z_selection_gal = (z_gal >= zmin) & (z_gal < zmax)
		N_gal = n.histogram(logm_gal[z_selection_gal], bins=bins_SMF)[0]

		p.figure(2, (6, 6))
		p.axes([0.16, 0.15, 0.8, 0.8])

		dx = n.log10(0.6777**2) 

		if z_mean[0]<=0.35:
			x_41, y_41, y_41_up, y_41_low = n.loadtxt( os.path.join( agn_data_dir, 'duty_cycle_G11_z01_LXhardgt41.ascii'), unpack=True)
			p.fill_between(x_41+dx, y1=10**(y_41_low),  y2=10**(y_41_up), label=r'G11 z=0.1 $L_X>10^{41}$erg s$^{-1}$', alpha=0.5, color='green')
			
		if z_mean[0]<=0.35:
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_S08_z025_LXhardgt42.ascii'), unpack=True)
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'S08 z=0.25 $L_X>10^{42}$erg s$^{-1}$', alpha=0.5, color='brown')
			
		if z_mean[0]<=0.5:
			x_41, y_41, y_41_up, y_41_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z025_LXhardgt41.ascii'), unpack=True)
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z025_LXhardgt42.ascii'), unpack=True)
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z025_LXhardgt43.ascii'), unpack=True)
			lab_bib = 'G17 z=0.25' 
			#p.fill_between(x_44, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r' $L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
			p.fill_between(x_41+dx, y1=10**(y_41_low),  y2=10**(y_41_up), label=r'$L_X>10^{41}$erg s$^{-1}$', alpha=0.5,color='black')
			
		if z_mean[0]>0.5 and z_mean[0]<1. :
			x_41, y_41, y_41_up, y_41_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt41.ascii'), unpack=True)
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt42.ascii'), unpack=True)
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt43.ascii'), unpack=True)
			x_44, y_44, y_44_up, y_44_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z075_LXhardgt44.ascii'), unpack=True)
			lab_bib = 'G17 z=0.75' 
			p.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
			p.fill_between(x_41+dx, y1=10**(y_41_low),  y2=10**(y_41_up), label=r'$L_X>10^{41}$erg s$^{-1}$', alpha=0.5,color='black')
			
		if z_mean[0]>1. and z_mean[0]<1.5 :
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z125_LXhardgt42.ascii'), unpack=True)
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z125_LXhardgt43.ascii'), unpack=True)
			x_44, y_44, y_44_up, y_44_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z125_LXhardgt44.ascii'), unpack=True)
			lab_bib = 'G17 z=1.25' 
			p.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
			#p.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')
			
		if z_mean[0]>1.5 and z_mean[0]<2. :
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z175_LXhardgt42.ascii'), unpack=True)
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z175_LXhardgt43.ascii'), unpack=True)
			x_44, y_44, y_44_up, y_44_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z175_LXhardgt44.ascii'), unpack=True)
			lab_bib = 'G17 z=1.75' 
			p.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
			#p.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')
			
		if z_mean[0]>2. and z_mean[0]<2.5 :
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z225_LXhardgt42.ascii'), unpack=True)
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z225_LXhardgt43.ascii'), unpack=True)
			x_44, y_44, y_44_up, y_44_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z225_LXhardgt44.ascii'), unpack=True)
			lab_bib = 'G17 z=2.25' 
			p.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
			#p.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')
			
		if z_mean[0]>2.5 and z_mean[0]<3. :
			x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z275_LXhardgt42.ascii'), unpack=True)
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z275_LXhardgt43.ascii'), unpack=True)
			x_44, y_44, y_44_up, y_44_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z275_LXhardgt44.ascii'), unpack=True)
			lab_bib = 'G17 z=2.75' 
			p.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			p.fill_between(x_42+dx, y1=10**(y_42_low),  y2=10**(y_42_up), label=r'$L_X>10^{42}$erg s$^{-1}$', alpha=0.5,color='red')
			#p.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')
			
		if z_mean[0]>3 :
			lab_bib = 'G17 z=3.5' 
			x_43, y_43, y_43_up, y_43_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z350_LXhardgt43.ascii'), unpack=True)
			x_44, y_44, y_44_up, y_44_low = n.loadtxt( os.path.join( agn_data_dir,  'duty_cycle_G17_z350_LXhardgt44.ascii'), unpack=True)
			#x_42, y_42, y_42_up, y_42_low = n.loadtxt( os.path.join( agn_data_dir, 'duty_cycle_H10_z02_LXhardgt42.ascii'), unpack=True)
			p.fill_between(x_44+dx, y1=10**(y_44_low),  y2=10**(y_44_up), label=lab_bib+r' $L_X>10^{44}$erg s$^{-1}$', alpha=0.5, color='magenta')
			p.fill_between(x_43+dx, y1=10**(y_43_low),  y2=10**(y_43_up), label=r'$L_X>10^{43}$erg s$^{-1}$', alpha=0.5,color='blue')
			#p.fill_between(x_42, y1=10**(y_42_low),  y2=10**(y_42_up),label=r'$L_X>10^{42}$', alpha=0.5,color='red')
			#p.fill_between(x_41, y1=10**(y_41_low),  y2=10**(y_41_up),label=r'$L_X>10^{41}$', alpha=0.5,color='black')

		p.axhline(f_duty(z_mean[0]), ls='dashed', lw=2)

		tsel = (z_selection)
		N_agn_a = n.histogram(logm[tsel], bins=bins_SMF)[0]
		y = N_agn_a * 1. / N_gal  # *dc_val
		yerr = y * N_agn_a**(-0.5)  # *dc_val
		p.errorbar(x_SMF, y, yerr=yerr, color='grey', label='all AGN')

		tsel = (lx > 41) & (z_selection)
		N_agn_41 = n.histogram(logm[tsel], bins=bins_SMF)[0]
		y = N_agn_41 * 1. / N_gal  # *dc_val
		yerr = y * N_agn_41**(-0.5)  # *dc_val
		p.errorbar(x_SMF, y, yerr=yerr, color='black', label=r'L$_X>10^{41}$')

		tsel = (lx > 42) & (z_selection)
		N_agn_42 = n.histogram(logm[tsel], bins=bins_SMF)[0]
		y = N_agn_42 * 1. / N_gal  # *dc_val
		yerr = y * N_agn_42**(-0.5)  # *dc_val
		p.errorbar(x_SMF, y, yerr=yerr, color='red', label=r'L$_X>10^{42}$')

		tsel = (lx > 43) & (z_selection)
		N_agn_43 = n.histogram(logm[tsel], bins=bins_SMF)[0]
		y = N_agn_43 * 1. / N_gal  # *dc_val
		yerr = y * N_agn_43**(-0.5)  # *dc_val
		p.errorbar(x_SMF, y, yerr=yerr, color='blue', label=r'L$_X>10^{43}$')

		tsel = (lx > 44) & (z_selection)
		N_agn_44 = n.histogram(logm[tsel], bins=bins_SMF)[0]
		y = N_agn_44 * 1. / N_gal  # *dc_val
		yerr = y * N_agn_44**(-0.5)  # *dc_val
		p.errorbar(x_SMF, y, yerr=yerr, color='magenta', label=r'L$_X>10^{44}$')

		p.xlabel(r'$\log_{10}(M^*/M_\odot)$')
		p.ylabel(r'$f_{AGN}(M^*, ' + str(n.round(zmin, 2)) + r"<z<" + str(n.round(zmax, 2)) + r')$')
		p.yscale('log')
		p.ylim((5e-5, 0.4))
		p.xlim((9.5, 12.))
		p.grid()
		# , '+str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)))
		p.title('Duty cycle')
		#p.legend( loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, title=str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)) )
		p.legend(frameon=False, loc=0)
		#p.legend( loc=8, ncol=2, fancybox=True, title=str(n.round(zmin,2))+"<z<"+str(n.round(zmax,2)) )
		p.savefig(os.path.join(DC_dir, "duty_cycle_AGN_" + baseName + ".png"))
		p.clf()

		p.figure(1, (6, 6))
		p.axes([0.17, 0.15, 0.73, 0.73])
		for fun, name in zip(smf_ilbert_fun, smf_ilbert_name):
			p.plot(mbins, fun(10**mbins)/(0.7**3), label=name, ls='dashed', lw=0.5)
		p.plot(x_SMF, N_gal / (vol * dlogM), color='green', label='all galaxies')
		p.plot(x_SMF, N_agn_a / (vol * dlogM), color='grey', label='all AGN')
		p.plot(x_SMF, N_agn_41 / (vol * dlogM), color='black', label=r'L$_X>10^{41}$')
		p.plot(x_SMF, N_agn_42 / (vol * dlogM), color='red', label=r'L$_X>10^{42}$')
		p.plot(x_SMF, N_agn_43 / (vol * dlogM), color='blue', label=r'L$_X>10^{43}$')
		p.plot(x_SMF, N_agn_44 / (vol * dlogM), color='magenta', label=r'L$_X>10^{44}$')
		p.xlabel(r'$\log_{10}(M^*/[M_\odot])$')
		p.ylabel(r'$\Phi$=dN/dlogL/dV [1/Mpc$^3$/dex]')
		p.legend(frameon=False, loc=0, fontsize=12)
		p.yscale('log')
		p.xlim((9, 12.5))
		p.ylim((1e-8, 1e-1))
		p.title('Stellar mass function, ' + str(n.round(zmin, 2)) + "<z<" + str(n.round(zmax, 2)))
		p.grid()
		p.savefig(os.path.join(DC_dir, "SMF_AGN_" + baseName + ".png"))
		p.clf()
