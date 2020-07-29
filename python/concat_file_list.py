"""
Concatenates a list of files

input: 
 - file containing a list of files to concatenate
 - path to the output file
example input : 

 ls $MD04/fits/cat_AGN_all_*/000000.fit > fit_list_MD04_all_000000_AGN_.list
 - p_2_file_list = "fit_list_MD04_all_000000_AGN_.list" 
 - path_2_output = "/data17s/darksim/simulation_3/MD/MD_0.4Gpc/cat_AGN_all/000000.fit" 
 
"""
import numpy as n
from astropy.table import Table, Column
import sys, os, glob 

p_2_file_list = sys.argv[1]
path_2_output = sys.argv[2]

file_list = n.loadtxt(p_2_file_list, dtype='str')

t = Table.read(file_list[0])
t['HALO_id'] = t['HALO_id'].astype('int64')     
t['HALO_pid'] = t['HALO_pid'].astype('int64')    
for p_2_file in file_list[1:]:
	print(p_2_file)
	t1 = Table.read(p_2_file)
	t1['HALO_id'] = t1['HALO_id'].astype('int64')     
	t1['HALO_pid'] = t1['HALO_pid'].astype('int64')    
	t = n.hstack(( t, t1))

t_out = Table(t)
t_out.write(path_2_output, overwrite=True)
