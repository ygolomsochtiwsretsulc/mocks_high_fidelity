"""

Pre-rpocess sixte event files to input eSASS

Adapted from bash to python

Based on the Script from Teng Liu

"""
import os, sys
import numpy as n

HID = sys.argv[1] # '767'                # sys.argv[1]
RA  = sys.argv[2] # '315.0'              #sys.argv[2]
DEC = sys.argv[3] # '-4.780191847199163' #sys.argv[3]
print( HID, RA, DEC )


#log_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass', 'logs-erass1')
log_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass', 'logs-erass8')
cluster_dir = os.path.join( "/data26s/comparat/simulations/erosim/eRASS8_cluster_MD10"   , HID.zfill(3) )
agn_dir   = os.path.join( "/data26s/comparat/simulations/erosim/eRASS8_agn_MD10"       , HID.zfill(3) )
stars_dir = os.path.join( "/data40s/erosim/eRASS/eRASS8_stars"                      , HID.zfill(3) )
pbg_dir   = os.path.join( "/data40s/erosim/eRASS/eRASS8_ParticleBackground", 'sim_particle_'+HID.zfill(3) )
bg_dir    = os.path.join( "/data40s/erosim/eRASS/diffuse_fg_rosat"                     , HID.zfill(3) )

pbg_file = os.path.join(pbg_dir, 'evt_particle_'+HID.zfill(3)+'.fits')

# where files are written :
#esass_dir   = os.path.join( "/data26s/comparat/simulations/erosim/SFC-v1-erass1" , HID.zfill(3) )
esass_dir   = os.path.join( "/data26s/comparat/simulations/erosim/SFC-v1-erass8" , HID.zfill(3) )

Attitude_File="/data40s/erosim/eRASS/eRASS_4yr_epc85_att.fits"

path_2_tmp_file = os.path.join(esass_dir, 'tmp_'+HID+'.fits')
path_2_event_file = os.path.join(esass_dir, 'evt_'+HID+'.fits')

if os.path.isdir(esass_dir)==False:
	os.mkdir( esass_dir )

def command_ero_cal_event( nccd = 1, topdir = agn_dir, prefix = 'a_', HID=HID, Attitude_File = Attitude_File ): 
	"""
	ero_calevents Projection=AIT Attitude=${ATTITUDE} clobber=yes EvtFile=${SIMU_DIR}/erass_ccd${Nccd}_evt.fits eroEvtFile=${SIMU_DIR}/cal_erass_ccd${Nccd}_evt.fits CCDNR=${Nccd} RA=${ra_cen} Dec=${dec_cen}
	"""
	part1 = "ero_calevents Projection=AIT Attitud="+Attitude_File
	part2 = " EvtFile="+os.path.join(topdir, 'erass_ccd' + str(nccd) + '_evt.fits') 
	part3 = " eroEvtFile="+os.path.join( esass_dir, prefix + HID +'_' + str(nccd) + '.fits') 
	part4 = " CCDNR="+str(nccd)+" RA="+RA+" DEC="+DEC+" clobber=yes"
	full_command = part1 + part2 + part3 + part4
	return full_command

def command_fparkey_15( nccd = 1, prefix = 'a_', HID=HID ): 
	"""
	fparkey 15 ${SIMU_DIR}/cal_erass_ccd${Nccd}_evt.fits[1] PAT_SEL add=yes
	"""
	return "fparkey 15 "+os.path.join( esass_dir, prefix + HID +'_' + str(nccd) + '.fits') +"[1] PAT_SEL add=yes"


def command_ftcopy( nccd = 1, prefix = 'a_', HID=HID ): 
	"""
	ftcopy ${SIMU_DIR}/cal_erass_ccd${Nccd}_evt.fits[1]"[col *,FRAMETIME=TIME,RECORDTIME=TIME]" ${SIMU_DIR}/evt_erass_ccd${Nccd}.fits clobber=yes

	"""
	part1 = """ftcopy '"""+os.path.join( esass_dir, prefix + HID +'_' + str(nccd) + '.fits')
	part2 = """[col *,RECORDTIME=TIME,FRAMETIME=TIME]'""" 
	part3 = " "+os.path.join( esass_dir, prefix + HID +'_' + str(nccd) + '_FTCOPY.fits')
	part4 = " clobber=yes"
	full_command = part1 + part2 + part3 + part4
	return full_command


if os.path.isfile(path_2_tmp_file)==False:
	for NCCD in n.arange(7)+1:
		print(NCCD)
		# a: AGN
		cmd_a = command_ero_cal_event(NCCD, agn_dir, 'a_')
		print(cmd_a)
		os.system(cmd_a)
		cmd_a_1 = command_fparkey_15(NCCD, 'a_')
		print(cmd_a_1)
		os.system(cmd_a_1)
		cmd_a_2 = command_ftcopy(NCCD, 'a_')
		print(cmd_a_2)
		os.system(cmd_a_2)
		# b: Background
		cmd_b = command_ero_cal_event(NCCD, bg_dir, 'b_')
		print(cmd_b)
		os.system(cmd_b)
		cmd_b_1 = command_fparkey_15(NCCD, 'b_')
		print(cmd_b_1)
		os.system(cmd_b_1)
		cmd_b_2 = command_ftcopy(NCCD, 'b_')
		print(cmd_b_2)
		os.system(cmd_b_2)
		# c: Clusters
		cmd_c = command_ero_cal_event(NCCD, cluster_dir, 'c_')
		print(cmd_c)
		os.system(cmd_c)
		cmd_c_1 = command_fparkey_15(NCCD, 'c_')
		print(cmd_c_1)
		os.system(cmd_c_1)
		cmd_c_2 = command_ftcopy(NCCD, 'c_')
		print(cmd_c_2)
		os.system(cmd_c_2)
		# d: stars
		cmd_d = command_ero_cal_event(NCCD, cluster_dir, 'd_')
		print(cmd_d)
		os.system(cmd_d)
		cmd_d_1 = command_fparkey_15(NCCD, 'd_')
		print(cmd_d_1)
		os.system(cmd_d_1)
		cmd_d_2 = command_ftcopy(NCCD, 'd_')
		print(cmd_d_2)
		os.system(cmd_d_2)



	path_2_list_input = os.path.join(log_dir, HID+'_input.list')
	cmd_list = "ls "+esass_dir+"/?_"+HID+"_[1-7].fits > "+path_2_list_input
	print( cmd_list )
	os.system( cmd_list )
	
	to_remove_input = n.loadtxt(path_2_list_input, dtype='str')
	print('number of files=',len(to_remove_input))
	for element in to_remove_input:
		os.remove(element)
	os.remove(path_2_list_input)
	
	# concatenates events 
	path_2_list = os.path.join(log_dir, HID+'_0.list')
	cmd_list = "ls "+esass_dir+"/?_"+HID+"_[1-7]_FTCOPY.fits "+pbg_file+" > "+path_2_list
	print( cmd_list )
	os.system( cmd_list )

	# evtool command (requires eSASS to be sourced)
	# evtool eventfiles="@tmp_${hid}.list" outfile=tmp_${hid}.fits pattern=15 clobber=yes
	cmd_evtool = """evtool eventfiles='@"""+path_2_list+"""' outfile="""+ path_2_tmp_file+" pattern=15 clobber=yes "
	#cmd_evtool_erass1 = """evtool eventfiles='@"""+path_2_list+"""' GTI="617943605 649479605" outfile="""+ path_2_tmp_file+" pattern=15 clobber=yes "
	print(cmd_evtool)
	os.system(cmd_evtool)

 
# gives the healpy center position
# radec2xy file=tmp_${hid}.fits ra0=${ra_cen} dec0=${dec_cen}
cmd_radec = "radec2xy file="+path_2_tmp_file+" ra0="+RA+" dec0="+DEC
print(cmd_radec)
os.system(cmd_radec)

# ftcopy tmp_${hid}.fits"[(RAWX-192.5)**2+(RAWY-192.5)**2<=191.5**2]" evt_${hid}.fits clobber=yes
cmd_copy = "ftcopy "+path_2_tmp_file+""""[(RAWX-192.5)**2+(RAWY-192.5)**2<=191.5**2]" """+path_2_event_file+" clobber=yes"
print(cmd_copy)
os.system(cmd_copy)

# fparkey SIXTE evt_${hid}.fits[1] CREATOR add=yes
cmd_parkey = "fparkey SIXTE "+path_2_event_file+"[1] CREATOR add=yes"
print(cmd_parkey)
os.system(cmd_parkey)

os.system('chgrp erosim '+esass_dir+'/*')
#to_remove = n.loadtxt(path_2_list, dtype='str')
#for element in to_remove:
	#os.remove(element)
#os.remove(path_2_list)
#os.remove(path_2_tmp_file)
