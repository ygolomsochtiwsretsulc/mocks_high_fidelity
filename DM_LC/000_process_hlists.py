"""

Replicates the snapshot over the sky and extracts halos in a given shell of comoving distance
Does full sky replication

Parameters
==========

 * l_box: length of the box in Mpc/h, example 400
 * h: unitless hubble parameter h=H0/100mk/s/Mpc, example 0.6777
 * N_skip: how many lines of header should be skipped while reading the halo file
 * env: environment variable linking to the simulation directory, example MD04 MD10, UNIT_fA1_DIR ...
 * path_2_header: path to the header of the output file
 * path_2_in_0: path to the rockstar hlist file, example /path/2/simulation/hlist_1.00000.list.bz2

Processing
==========

bunzip the file (optional)

Generates and executes awk commands to replicate and filter the halo files. Lines 120 - 129 are the most important lines of the script.

Concatenate all replicated files with a header

Convert to a fits file

outputs
=======

 * directory(path_2_in_0)/replication_$expansionParameter/all the files: contains a single file per replication. Deleted after the computation finishes.
 * directory(path_2_in_0)/replicated_$expansionParameter/one summary file: contains one file with all replication concatenated. Could be deleted after the computation is done.
 * $env/fits/all_'+$expansionParameter+'.fits') the former file converted to fits staged into the final directory. Keep this one ! It is the shell containing all the halos

"""
print('runs 001_process_hlists.py with arguments' )
import time
t0=time.time()
import numpy as n
import os, sys, subprocess
print(sys.argv)

l_box = float(sys.argv[1]) # = 1000.0 
h_str = sys.argv[2] # = 0.6777
N_skip = int(sys.argv[3])
env = sys.argv[4]
Mmin_str = sys.argv[5] # '9.0e9'
path_2_in_0 = sys.argv[6]

# deduce the L_box without h
h = float(h_str)
L_box = l_box/h
path_2_header =os.path.join( os.environ[env], 'header')

# if the file is zipped, first unzip
# else does nothing
# unzipped name is assigned in the variable: path_2_in
if path_2_in_0.split('.')[-1] == 'bz2':
    print('bunzip2 the file')
    os.chdir(os.path.dirname(path_2_in_0))
    os.system('bunzip2 ' + path_2_in_0)
    path_2_in = path_2_in_0[:-4]
else:
    path_2_in = path_2_in_0

# reads the list of snapshot available and the boundaries to be applied for each
# previously computed with 000_geometry.py
N_snap, Z_snap, A_snap, DC_max, DC_min, NN_replicas = n.loadtxt(os.path.join(
    os.environ[env], 'snap_list_with_border.txt'), unpack=True)


def get_a(path_2_in):
    """
    Function to retrieves the expansion parameter from the name of the file, not form the header as headers can vary from version to version.
    """
    alp = os.path.basename(path_2_in).split('_')[1][:7]
    print('a=', alp)
    return float(alp)
    
# Retrieve the snapshot and gets the boundaries
a_snap = get_a(path_2_in)
sel = (A_snap == a_snap)
# Get the NN parameter: very important
# NN: number of times the box is replicated by 1 Lbox positively and negatively. Example: NN=3 , the box gets replicated by -3xL_box, -2xL_box, -1xL_box, 0xL_box, 1xL_box, 2xL_box in all X,Y,Z directions, (2*NN)^3 = 6^3 = 216 replications.
NN = NN_replicas[sel][0]

# creates the directory where the outputs will be temporarily stored :
out_dir = os.path.join(
	os.path.dirname(path_2_in),
	'replication_' + str(a_snap) + '/')
if os.path.isdir(out_dir) == False:
	os.mkdir(out_dir)

out_dir_2 = os.path.join(
	os.path.dirname(path_2_in),
	'replicated_' + str(a_snap) + '/')
if os.path.isdir(out_dir_2) == False:
	os.mkdir(out_dir_2)

# Creates the regular pavement for the replication
pts_i = n.arange(-1 * NN, NN, 1)
ix_i, iy_i, iz_i = n.meshgrid(pts_i, pts_i, pts_i)
ix, iy, iz = n.ravel(ix_i), n.ravel(iy_i), n.ravel(iz_i)
# print(ix,iy,iz)
pts = pts_i * L_box
x0_i, y0_i, z0_i = n.meshgrid(pts, pts, pts)
x0, y0, z0 = n.ravel(x0_i), n.ravel(y0_i), n.ravel(z0_i)

######################
# WRITES and executes THE AWK COMMANDS
######################
# N_skip=64
def transform_rockstar_out_catalog_into_small_ascii_catalog_distinct(path_2_in, path_2_out, xt, D2_min, D2_max, Mmin_str=Mmin_str, N_skip=N_skip):
    """
    Replicates one time the simulation with offset x0,y0,z0=xt
    Keeps halos with D2_min^2<x^2+y^2+z^2<D2_max^2
    
    writes out only the interesting columns and the translated x,y,z coordinates for elements present in the shell DC_min, DC_max
    
    parameters
    
    * path_2_in: path to input files
    * path_2_out: path to output files
    * xt: offset vector, contains x0,y0,z0
    * D2_min: minimum distance
    * D2_max: maximum distance
    * Mmin_str: minimum virial mass string
    * N_skip: number of lines skipped at the beginning of the file (header), typically 64 for the Rockstar halo files
    
    """
    x0, y0, z0 = xt
    gawk_start = """gawk 'NR>""" + str(N_skip)
    # mvir mass filter
    gawk_filter_mass = """ {if ( ( $11 >= """+Mmin_str+""" ) """
    # minimum distance (shell)
    gawk_filter_dmin = """&& ( ($18/"""+h_str+"""+""" + str(x0) + """) ** 2 + ($19/"""+h_str+"""+""" + str(y0) + """) ** 2 + ($20/"""+h_str+"""+""" + str(z0) + """) ** 2 >=  """ + str(D2_min) + """ ) """
    # maximum distance (shell)
    gawk_filter_dmax = """&& ( ($18/"""+h_str+"""+""" + str(x0) + """) ** 2 + ($19/"""+h_str+"""+""" + str(y0) + """) ** 2 + ($20/"""+h_str+"""+""" + str(z0) + """) ** 2 <  """ + str(D2_max) + """ ) """
    # print columns
    gawk_print = """) print  $2, $6, $11, $12, $13, $16, $17, $18/"""+h_str+"""+""" + str(x0) + """, $19/"""+h_str+"""+""" + str(y0) + """, $20/"""+h_str+"""+""" + str(z0) + """, $21, $22, $23, $41, $42, $44, $52, $53, $67}' """
    gawk_end = path_2_in + " > " + path_2_out
    # full command :
    gawk_command = gawk_start + gawk_filter_mass + gawk_filter_dmax + gawk_filter_dmin + gawk_print + gawk_end
    print(gawk_command)
    print('-------------------------------')
    os.system(gawk_command)

# Loop over all (2*NN)^3 configuration to do all the replications
XT = n.transpose([x0, y0, z0])
iXT = n.transpose([ix, iy, iz])
for jj in range(len(iXT)):
    jx,jy,jz=iXT[jj]
    path_2_out = os.path.join(  out_dir, 'replication_'+str(jx)+'_'+str(jy)+'_'+str(jz)+'_hlist_'+str(a_snap)+'.list')
    transform_rockstar_out_catalog_into_small_ascii_catalog_distinct(path_2_in, path_2_out, XT[jj], DC_min[sel][0]**2, DC_max[sel][0]**2, Mmin_str=Mmin_str, N_skip=N_skip)

# list all replicated boxes 
# compute file size
# delete the ones that are empty
out = os.listdir(out_dir)
fs = n.array([ os.path.getsize(os.path.join( out_dir,oo)) for oo in out ])
to_delete = n.array(out)[(fs==0)]
n.array([ os.remove(os.path.join(out_dir, tdl)) for tdl in to_delete ])

# define location of the concatenated file
out_file = os.path.join( out_dir_2, 'all_'+str(a_snap)+'.list')
# list again all replicated boxes (the empty ones were deleted)
all_files = n.array([ os.path.join(out_dir, el) for el in n.array(out)[(fs>0)] ])

# bash command cat to concatenate the files with the header 
to_concat = " ".join(all_files)
command = 'cat '+path_2_header+' ' + to_concat + ' > '+ out_file
os.system(command)

# define location of the concatenated fits file
fits_dir = os.path.join( os.environ[env], 'fits')
if os.path.isdir(fits_dir) == False:
    os.mkdir(fits_dir)

out_fits = os.path.join(fits_dir, 'all_' + os.path.basename(path_2_in).split('_')[-1])
# convert the ascii file to the fits format 
stilts_command = "stilts tpipe in="+out_file+" ifmt=ascii omode=out ofmt=fits out="+out_fits[:-5]+".fits"
print(stilts_command)
os.system(stilts_command)

print('deletes intermediate files')
os.system('rm -rf '+ out_dir_2)
os.system('rm -rf '+ out_dir)

print('bzip2 ', path_2_in)
os.chdir(os.path.dirname(path_2_in))
os.system('bzip2 ' + path_2_in)

print('finished, dt=', time.time()-t0, 'seconds')
