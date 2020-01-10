"""
Writes a set of spectra in the simput format
"""
import xspec
import numpy as n
import os, sys
import astropy.io.fits as fits
xspec.Xset.cosmo = "67.77 0. 0.692885"

#s1 = xspec.Spectrum("arf01_100nmAl_200nmPI_sdtq.fits")
#b1 = s1.background
#r1 = s1.response
#arfFileName = r1.arf

#torus_dir = os.path.join(os.environ['DARKSIM_DIR'], 'model', 'torus_buchner' )
#torus_model = os.path.join(torus_dir, 'uxclumpy-cutoff.fits')
#torus_model_omni = os.path.join(torus_dir, 'uxclumpy-cutoff-omni.fits')
#torus_model_reflect = os.path.join(torus_dir, 'uxclumpy-cutoff-reflect.fits')
#torus_model_transmit = os.path.join(torus_dir, 'uxclumpy-cutoff-transmit.fits')


nh_val = 1
PL=1.9 
redshift=0.
f_scatter=0.02
norm1 = 1-f_scatter
norm2 = f_scatter
norm3 = 1.
rel_refl= -1.
incl = n.cos(30.*n.pi/180.)
norm_PR = 1.


def create_fits_file_spectrum(out_name, nH, redshift, energies, flux_density):
	n_values = str(len(energies))
	print(energies.shape, flux_density.shape)
	print(n.array([energies], dtype='float'))
	col1 = fits.Column(name='ENERGY', unit='keV', format = '1PE('+n_values+')', array=[energies] )
	col2 = fits.Column(name='FLUXDENSITY', unit='photon/s/cm**2/keV', format = '1PE('+n_values+')', array=[flux_density] )
	cols = fits.ColDefs([col1, col2])
	#tbhdr = fits.Header()
	#tbhdr['EXTNAME'] = 'SPECTRUM'
	#tbhdr['HDUCLASS'] = 'HEASARC/SIMPUT'
	#tbhdr['HDUCLAS1'] = 'SPECTRUM'
	#tbhdr['HDUVERS'] = '1.1.0'
	#tbhdr['EXTVER'] = 1 # or 0
	hdu = fits.BinTableHDU.from_columns(cols)
	hdu.name = 'SPECTRUM'
	hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
	hdu.header['HDUCLAS1'] = 'SPECTRUM'
	hdu.header['HDUVERS'] = '1.1.0'
	hdu.header['RADESYS'] = 'FK5'
	hdu.header['EQUINOX'] = 2000.0
	hdu.header['nH'] = nH # or 0
	hdu.header['z'] = redshift # or 0
	
	outfile = os.path.join(os.environ['MD10'], 'sixte', 'agn', 'spectra', out_name+'.fits')
	if os.path.isfile(outfile):
		os.system("rm "+outfile)
	hdu.writeto(outfile , clobber=True)

def get_spectrum(nh_val, redshift, d_kev = 0.0025):
	m1 = xspec.Model("(tbabs*(plcabs+pexrav)+zpowerlw)*tbabs")
	m1.pexrav.rel_refl='-2 -2 -2 -2'
	m1.setPars(
		nh_val,   #1    1   TBabs      nH         10^22    1.00000      +/-  0.0          
		nh_val,   #2    2   plcabs     nH         10^22    1.00000      +/-  0.0          
		3.,       #3    2   plcabs     nmax       (scale)  1.00000      
		1.,       #4    2   plcabs     FeAbun              1.00000      frozen
		7.11,     #5    2   plcabs     FeKedge    KeV      7.11000      frozen
		PL,       #6    2   plcabs     PhoIndex            2.00000      +/-  0.0          
		50.,      #7    2   plcabs     HighECut   keV      95.0000      frozen
		200.,     #8    2   plcabs     foldE               100.000      frozen
		1.0,      #9    2   plcabs     acrit               1.00000      frozen
		0.0,      #10    2   plcabs     FAST       (scale)  0.0          
		redshift, #11    2   plcabs     Redshift            0.0          frozen
		norm1,    #12    2   plcabs     norm                1.00000      +/-  0.0          
		PL,       #13    3   pexrav     PhoIndex            2.00000      +/-  0.0          
		200.,     #14    3   pexrav     foldE      keV      100.000      +/-  0.0          
		rel_refl, #15    3   pexrav     rel_refl            0.0          +/-  0.0          
		redshift, #16    3   pexrav     Redshift            0.0          frozen
		1.,       #17    3   pexrav     abund               1.00000      frozen
		1.,       #18    3   pexrav     Fe_abund            1.00000      frozen
		incl,     #19    3   pexrav     cosIncl             0.450000     frozen
		norm_PR,  #20    3   pexrav     norm                1.00000      +/-  0.0          
		PL,       #21    4   zpowerlw   PhoIndex            1.00000      +/-  0.0          
		redshift, #22    4   zpowerlw   Redshift            0.0          frozen
		norm2, #23    4   zpowerlw   norm                1.00000      +/-  0.0          
		0.01 #24    5   TBabs      nH         10^22    1.00000      +/-  0.0          
		)
	m1.pexrav.norm.link='p12/(1. + p11)/(1./(1. + p11))^( - p6)'
	#kevs = n.arange(0.1, 50, 0.1) #
	kevs = 10**n.arange(-1.,n.log10(50),d_kev)
	fluxes = []
	nPh = []
	for kev_min_erosita_RF, kev_max_erosita_RF in zip(kevs[:-1],kevs[1:]):
		xspec.AllModels.calcFlux(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
		#print(xspec.Xset.cosmo)
		#xspec.AllModels.calcLumin(str(kev_min_erosita_RF)+" "+str(kev_max_erosita_RF))
		fluxes.append( m1.flux[0]    )
		nPh.append(m1.flux[3] )
	return (kevs[:-1]+kevs[1:])*0.5, n.array(fluxes), n.array(nPh), -kevs[:-1]+kevs[1:]

def create_spectrum(redshift, nh_val, d_kev = 0.0025):
	energies,y,nP,dE=get_spectrum(nh_val, redshift, d_kev = d_kev)
	n_values = str(len(energies))
	out_name = 'agn_nH_'+str(n.round(n.log10(nh_val*10**22),1))+'_z_'+str(n.round(redshift,1))+'_nEbins_'+n_values
	print(out_name)
	flux_density = nP/dE
	create_fits_file_spectrum(out_name, nh_val, redshift, energies, flux_density)


nh_vals = 10**n.arange(-2,4+0.01,0.2)#0.05)
z_vals = n.arange(0., 6.01, 0.1)
KEVS = n.array([0.32, 0.16, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025])
#KEVS = n.array([0.64])

NH,ZZ = n.meshgrid(nh_vals, z_vals)
for d_kev in KEVS:
	for nh_val, z_val in zip(n.hstack((NH)),n.hstack((ZZ))) :
		create_spectrum( z_val, nh_val, d_kev )




