"""
What it does
------------

Creates a simput catalog fo each healpix pixel (healpix pixel of 13.4 deg2)

References
----------

Command to run
--------------

python3 004_6_clusters_simput.py environmentVAR

arguments
---------

environmentVAR: environment variable linking to the directory where files are e.g. "MD10"
It will then work in the directory : $environmentVAR/hlists/fits/

Dependencies
------------

import time, os, sys, glob, numpy, astropy, scipy, matplotlib

/data40s/erosim/eRASS/eRASS8_cluster/000/erass_ccd1_evt.fits
/data17s/darksim/MD/MD_1.0Gpc/cat_CLU_SIMPUT/
"""
from astropy.table import Table, Column
import sys, os, time
import astropy.units as u
from astropy_healpix import healpy
import astropy.io.fits as fits
import numpy as n
print('CREATES SIMPUT CLUSTER FILES')
print('------------------------------------------------')
print('------------------------------------------------')
t0 = time.time()

#from sklearn.neighbors import BallTree

#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 14})
#import matplotlib.pyplot as p
print(sys.argv)
env = sys.argv[1] #
pixel_size_image = float(sys.argv[2]) # 20 or 2
b_HS = float(sys.argv[3]) 
cov_mat_option = sys.argv[4] 

# initializes pathes to files

root_dir = os.path.join(os.environ[env])

path_2_catalog = os.path.join(root_dir, env+'_eRO_CLU_b'+str(int(b_HS*10))+'_CM_'+cov_mat_option + "_pixS_" + str(pixel_size_image)+'.fits')
dir_2_eRO_all = os.path.join(root_dir, 'cat_eRO_CLU_b'+str(int(b_HS*10))+'_CM_'+cov_mat_option + "_pixS_" + str(pixel_size_image))
dir_2_SMPT = os.path.join(root_dir, 'cat_CLU_SIMPUT_b'+str(int(b_HS*10))+'_CM_'+cov_mat_option + "_pixS_" + str(pixel_size_image))
dir_2_SMPT_image = os.path.join(dir_2_SMPT, 'cluster_images')
dir_2_SMPT_spectra = os.path.join(dir_2_SMPT, 'cluster_Xspectra')

print(path_2_catalog)
print(dir_2_eRO_all)
print(dir_2_SMPT)
print(dir_2_SMPT_image)
print(dir_2_SMPT_spectra)

if os.path.isdir(dir_2_eRO_all) == False:
    os.system('mkdir -p ' + dir_2_eRO_all)
if os.path.isdir(dir_2_SMPT_image) == False:
    os.system('mkdir -p ' + dir_2_SMPT_image)
if os.path.isdir(dir_2_SMPT_spectra) == False:
    os.system('mkdir -p ' + dir_2_SMPT_spectra)
if os.path.isdir(dir_2_SMPT) == False:
    os.system('mkdir -p ' + dir_2_SMPT)

t = Table.read(path_2_catalog)
HEALPIX_8 = healpy.ang2pix(8, n.pi/2. - t['DEC']*n.pi/180. , t['RA']*n.pi/180. , nest=True)

"""
# subtract flux emmitted by satellites from the distinct halos that host them 
out = n.unique(t['HALO_id'].astype('int'), return_counts=True)
out2 = n.unique(t['HALO_pid'].astype('int'), return_counts=True)
cen = (t['HALO_pid']==-1)
sat = (t['HALO_pid']>=0) & (t['CLUSTER_FX_soft'] > 10**(t['flux_limit_eRASS8_pt']))

out_sat = n.unique(t['HALO_pid'][sat], return_counts=True, return_index=True, return_inverse=True) 
out_cen = n.unique(t['HALO_id'][cen], return_counts=True, return_index=True, return_inverse=True) 
host_ids = out_sat[0]
for el in host_ids:
	subss = (t['HALO_pid'] == el )
	hosts = (t['HALO_id']  == el )
	to_remove = t['CLUSTER_FX_soft'][subss]
	t_host = t['CLUSTER_FX_soft'][hosts]
	if len(to_remove)==1 and len(t_host)==len(to_remove):
		val_all = t['CLUSTER_FX_soft'][hosts]
		val_new = val_all - to_remove 
		t['CLUSTER_FX_soft'][hosts] = val_new
	elif len(to_remove)>1 and len(t_host)==len(to_remove):
		#print('Nhost',len(t_host),',Nsat',len(to_remove), 'id=', el)
		Tree_profiles = BallTree(n.transpose([t['RA'][hosts], t['DEC'][hosts]]))
		DATA = n.transpose([t['RA'][subss], t['DEC'][subss]])
		ids_out = Tree_profiles.query(DATA, k=1, return_distance = False)
		ids = n.hstack((ids_out))
		val_all = t['CLUSTER_FX_soft'][hosts][ids]
		val_new = val_all - to_remove 
		t_host = t['CLUSTER_FX_soft'][hosts]
		t_host[ids] = val_new
		t['CLUSTER_FX_soft'][hosts] = t_host
	elif len(to_remove)>1 and len(t_host)<len(to_remove) and len(t_host)>=1:
		#print('SingleHost, Nhost',len(t_host),',Nsat',len(to_remove), 'id=', el)
		val_all = t['CLUSTER_FX_soft'][hosts]
		val_new = val_all - to_remove.sum() 
		t['CLUSTER_FX_soft'][hosts] = val_new
	else :
		print('Orphans, Nhost',len(t_host),',Nsat',len(to_remove), 'id=', el)
""" 


# link to templates
#def tpl_name(temperature, redshift): return 'cluster_Xspectra/cluster_spectrum_10kT_' + str(int(temperature * 1000)).zfill(4) + '_100z_' + str(int(redshift * 10000)).zfill(4) + '.fits[SPECTRUM][#row==1]'
def tpl_name(temperature, redshift): return 'cluster_Xspectra/cluster_spectrum_10000kT_' + str(int(temperature * 10000)).zfill(7) + '_10000z_' + str(int(redshift * 10000)).zfill(7) + '.fits[SPECTRUM][#row==1]'
#HEALPIX_8_id = 151


for HEALPIX_8_id in n.arange(healpy.nside2npix(8)):
	"""
	Loops over healpix pixels and writes the files to path_2_eRO_catalog
	"""
	print(HEALPIX_8_id)
	path_2_SMPT_catalog = os.path.join(dir_2_SMPT, 'c_'+str(HEALPIX_8_id).zfill(6) + '.fit')
	# print(path_2_eRO_catalog)
	sf = (HEALPIX_8 == HEALPIX_8_id) # & ( t['redshift_R'] < 1.5 )
	t2 = t[sf]
	#if pixel_size_image>18:
		#bright = (t1['HALO_M500c'] > 7e13) # ( t1['CLUSTER_FX_soft'] > 10**( t1['flux_limit_eRASS8_pt'] - 0.3 ) ) # to get 300,000 over the full sky
	#else:
		#bright = ( t1['CLUSTER_FX_soft'] > 10**( -15 )  ) & ( HEALPIX_8[sf] == 0 ) # to get 300,000 over the full sky
	#t2 = t1[bright]
	N_clu_all = len(t2['RA'])
	print('Number of clusters=',N_clu_all)
	print('density=',N_clu_all/53.)
	if N_clu_all>0:
		ra_array = t2['RA']
		dec_array = t2['DEC']
		redshift = t2['redshift_R']
		FX_soft_attenuated = t2['CLUSTER_FX_soft']
		kT = t2['CLUSTER_kT']
		# size of the pixel in the image written
		# randomize orientations
		rd_all = n.random.rand(N_clu_all)
		orientation = n.random.rand(N_clu_all) * 180.  # IMGROTA
		# scale the image with the size of the cluster
		# all images have 5.5e-04*120*60 = 3.96 arc minute on the side
		# default size 0.033 arcmin/pixel
		pixel_rescaling =  n.ones_like(t2['angularSize_per_pixel']) 
		# NOW ASSIGNS TEMPLATES BASED ON THE HALO PROPERTIES
		#template = n.zeros(N_clu_all).astype('U100')
		#template[template == "0.0"] = "cluster_images/elliptical_ba_0p25_cc.fits[SPECTRUM][#row==1]"
		template = n.array([ el.strip()+".fits[IMAGE]" for el in t2['XRAY_image_path'] ])
		template_fullPath = n.array([ os.path.join( dir_2_SMPT, el.strip()+".fits")  for el in t2['XRAY_image_path'] ])
		template_exists = n.array([ os.path.isfile( os.path.join( dir_2_SMPT, el.strip()+".fits") ) for el in t2['XRAY_image_path'] ])
		N_exist = len(template_exists.nonzero()[0])
		N_templates = len(template_exists)
		if N_exist<N_templates:
			print('error, missing images', N_exist, N_templates)
		# NOW links to the grid of SPECTRA
		# kt_arr = n.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
		# z_arr = n.hstack((n.array([0., 0.05]), n.arange(0.1, 4., 0.1)))
		kt_arr = 10**n.arange(-1,1.3,0.01)
		z_arr = n.hstack((n.array([0., 0.01]), 10**n.arange(n.log10(0.02), n.log10(4.), 0.01)))
		# 
		indexes_kt = n.array([(n.abs(kT_val - kt_arr)).argmin() for kT_val in kT])
		kT_values = kt_arr[indexes_kt]
		indexes_z = n.array([(n.abs(z_val - z_arr)).argmin() for z_val in redshift])
		z_values = z_arr[indexes_z]
		spec_names = n.zeros(N_clu_all).astype('U200')
		# "cluster_Xspectra/cluster_spectrum_10kT_0100_100z_0150.fits[SPECTRUM][#row==1]"
		for jj, (kT_values_ii, z_values_ii) in enumerate(zip(kT_values, z_values)):
			spec_names[jj] = tpl_name(kT_values_ii, z_values_ii)

		N_per_simput = 999
		for jj, (id_min, id_max) in enumerate(zip(n.arange(0,N_clu_all,N_per_simput), n.arange(0,N_clu_all,N_per_simput)+N_per_simput)):
			path_2_SMPT_catalog = os.path.join(dir_2_SMPT, 'c_'+str(HEALPIX_8_id).zfill(6) + '_N_'+str(jj)+'.fit')
			hdu_cols = fits.ColDefs([
				fits.Column(name="SRC_ID",  format='K',    unit='',    array=(n.arange(N_clu_all) + 4e8).astype('int')[id_min:id_max]), 
				fits.Column(name="RA",      format='D',    unit='deg', array=ra_array[id_min:id_max]), 
				fits.Column(name="DEC",     format='D',    unit='deg', array=dec_array[id_min:id_max]), 
				fits.Column(name="E_MIN",   format='D',    unit='keV', array=n.ones(N_clu_all)[id_min:id_max] * 0.5), 
				fits.Column(name="E_MAX",   format='D',    unit='keV', array=n.ones(N_clu_all)[id_min:id_max] * 2.0), 
				fits.Column(name="FLUX",    format='D',    unit='erg/s/cm**2', array=FX_soft_attenuated[id_min:id_max]), 
				fits.Column(name="IMAGE",   format='100A', unit='', array=template[id_min:id_max]), 
				fits.Column(name="SPECTRUM",format='100A', unit='', array=spec_names[id_min:id_max]),         
				fits.Column(name="IMGROTA", format='D',    unit='deg', array=orientation[id_min:id_max]), 
				fits.Column(name="IMGSCAL", format='D',    unit='', array=pixel_rescaling[id_min:id_max])
			])
			hdu = fits.BinTableHDU.from_columns(hdu_cols)
			hdu.name = 'SRC_CAT'
			hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
			hdu.header['HDUCLAS1'] = 'SRC_CAT'
			hdu.header['HDUVERS'] = '1.1.0'
			hdu.header['RADESYS'] = 'FK5'
			hdu.header['EQUINOX'] = 2000.0
			outf = fits.HDUList([fits.PrimaryHDU(), hdu])  # ,  ])
			if os.path.isfile(path_2_SMPT_catalog):
				os.system("rm " + path_2_SMPT_catalog)
			outf.writeto(path_2_SMPT_catalog, overwrite=True)
			#print(path_2_SMPT_catalog, 'written', time.time() - t0)
			path_2_CLU_catalog = os.path.join(dir_2_eRO_all, 'c_'+str(HEALPIX_8_id).zfill(6) +'_N_'+str(jj)+ '.fit')
			t_out = Table( t2[id_min:id_max] )
			t_out.add_column(Column(name="SRC_ID",   dtype = n.int64,    unit='',            data = (n.arange(N_clu_all) + 4e8).astype('int')[id_min:id_max])   )
			t_out.add_column(Column(name="E_MIN",    dtype = n.float,    unit='keV',         data = n.ones(N_clu_all)[id_min:id_max] * 0.5)     )
			t_out.add_column(Column(name="E_MAX",    dtype = n.float,    unit='keV',         data = n.ones(N_clu_all)[id_min:id_max] * 2.0)     )
			t_out.add_column(Column(name="FLUX",     dtype = n.float,    unit='erg/s/cm**2', data = FX_soft_attenuated[id_min:id_max])          )
			t_out.add_column(Column(name="IMAGE",    dtype = n.str,      unit='',            data = template[id_min:id_max])                    )
			t_out.add_column(Column(name="SPECTRUM", dtype = n.str,      unit='',            data = spec_names[id_min:id_max])                  )
			t_out.add_column(Column(name="IMGROTA",  dtype = n.float,    unit='deg',         data = orientation[id_min:id_max])                 )
			t_out.add_column(Column(name="IMGSCAL",  dtype = n.float,    unit='',            data = pixel_rescaling[id_min:id_max])             )
			t_out.write(path_2_CLU_catalog, overwrite=True)
			#print(path_2_CLU_catalog, 'written', time.time() - t0)
