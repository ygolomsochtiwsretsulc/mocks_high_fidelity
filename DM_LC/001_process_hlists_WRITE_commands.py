#!/bin/bash
import numpy as n
import glob, os, sys

# UNIT_fA1_DIR

l_box = 1000.0 # float(sys.argv[1]) # = 1000.0 
h_str = 0.6774 # sys.argv[2] # = 0.6777
N_skip = 64 # int(sys.argv[3])
env = 'UNIT_fA1_DIR' # sys.argv[4]
Mmin_str = 2e11 # sys.argv[5] # '9.0e9'
#path_2_in_0 = sys.argv[6]

xx = sorted( n.array( glob.glob(os.path.join(os.environ[env],'hlists','hlist_?.?*'))))
params = str(l_box)+' '+str(h_str)+' '+str(N_skip)+' '+env+' '+str(Mmin_str)+' '
for x in xx[::-1]:
    cmd = "python 000_process_hlists.py " + params + x
    print(cmd)


# UNIT_fA1i_DIR

l_box = 1000.0 # float(sys.argv[1]) # = 1000.0 
h_str = 0.6774 # sys.argv[2] # = 0.6777
N_skip = 64 # int(sys.argv[3])
env = 'UNIT_fA1i_DIR' # sys.argv[4]
Mmin_str = 2e11 # sys.argv[5] # '9.0e9'
#path_2_in_0 = sys.argv[6]

xx = sorted( n.array( glob.glob(os.path.join(os.environ[env],'hlists','hlist_?.?*'))))
params = str(l_box)+' '+str(h_str)+' '+str(N_skip)+' '+env+' '+str(Mmin_str)+' '
for x in xx[::-1]:
    cmd = "python 000_process_hlists.py " + params + x
    print(cmd)


# UNIT_fA2_DIR

l_box = 1000.0 # float(sys.argv[1]) # = 1000.0 
h_str = 0.6774 # sys.argv[2] # = 0.6777
N_skip = 64 # int(sys.argv[3])
env = 'UNIT_fA2_DIR' # sys.argv[4]
Mmin_str = 2e11 # sys.argv[5] # '9.0e9'
#path_2_in_0 = sys.argv[6]

xx = sorted( n.array( glob.glob(os.path.join(os.environ[env],'hlists','hlist_?.?*'))))
params = str(l_box)+' '+str(h_str)+' '+str(N_skip)+' '+env+' '+str(Mmin_str)+' '
for x in xx[::-1]:
    cmd = "python 000_process_hlists.py " + params + x
    print(cmd)

# MD04

l_box = 400.0 # float(sys.argv[1]) # = 1000.0 
h_str = 0.6777 # sys.argv[2] # = 0.6777
N_skip = 64 # int(sys.argv[3])
env = 'MD04' # sys.argv[4]
Mmin_str = 1e10 # sys.argv[5] # '9.0e9'
#path_2_in_0 = sys.argv[6]

xx = sorted( n.array( glob.glob(os.path.join(os.environ[env],'hlists','hlist_?.?*'))))
params = str(l_box)+' '+str(h_str)+' '+str(N_skip)+' '+env+' '+str(Mmin_str)+' '
for x in xx[::-1]:
    cmd = "python 000_process_hlists.py " + params + x
    print(cmd)

# MD10

l_box = 1000.0 # float(sys.argv[1]) # = 1000.0 
h_str = 0.6777 # sys.argv[2] # = 0.6777
N_skip = 64 # int(sys.argv[3])
env = 'MD10' # sys.argv[4]
Mmin_str = 2e11 # sys.argv[5] # '9.0e9'
#path_2_in_0 = sys.argv[6]

xx = sorted( n.array( glob.glob(os.path.join(os.environ[env],'hlists','hlist_?.?*'))))
params = str(l_box)+' '+str(h_str)+' '+str(N_skip)+' '+env+' '+str(Mmin_str)+' '
for x in xx[::-1]:
    cmd = "python 000_process_hlists.py " + params + x
    print(cmd)

# MD40

l_box = 4000.0 # float(sys.argv[1]) # = 1000.0 
h_str = 0.6777 # sys.argv[2] # = 0.6777
N_skip = 64 # int(sys.argv[3])
env = 'MD40' # sys.argv[4]
Mmin_str = 1e13 # sys.argv[5] # '9.0e9'
#path_2_in_0 = sys.argv[6]

xx = sorted( n.array( glob.glob(os.path.join(os.environ[env],'hlists','hlist_?.?*'))))
params = str(l_box)+' '+str(h_str)+' '+str(N_skip)+' '+env+' '+str(Mmin_str)+' '
for x in xx[::-1]:
    cmd = "python 000_process_hlists.py " + params + x
    print(cmd)
