import os, sys, glob
import numpy as n

all_fields = n.arange(768)
# determine all_fields not in all_fields
index = []
for nccd in n.arange(7)+1:
	finished_i = n.array(glob.glob( "/data40s/erosim/eRASS/bkg_rosat/???/erass_ccd"+str(nccd)+"_evt.fits" ))
	finished = n.array([ int(el.split('/')[-2]) for el in finished_i ])
	index.append( n.in1d(all_fields, finished) )

index = n.array((index))
sum_index = n.sum(index, axis=0)
probably_completed_bkg = (sum_index==7)

all_fields = n.arange(768)
# determine all_fields not in all_fields
index = []
for nccd in n.arange(7)+1:
	finished_i = n.array(glob.glob( "/data40s/erosim/eRASS/bkg_rosat/???/erass_ccd"+str(nccd)+"_raw.fits" ))
	finished = n.array([ int(el.split('/')[-2]) for el in finished_i ])
	index.append( n.in1d(all_fields, finished) )

index = n.array((index))
sum_index = n.sum(index, axis=0)
not_completed_bkg = (sum_index==7)

completed_bkg = (probably_completed_bkg) & (not_completed_bkg==False)

index = []
for nccd in n.arange(7)+1:
	finished_i = n.array(glob.glob( "/data40s/erosim/eRASS/eRASS8_agn_MD10/???/erass_ccd"+str(nccd)+"_evt.fits" ))
	finished = n.array([ int(el.split('/')[-2]) for el in finished_i ])
	index.append( n.in1d(all_fields, finished) )

index = n.array((index))
sum_index = n.sum(index, axis=0)
completed_agn = (sum_index==7)

index = []
for nccd in n.arange(7)+1:
	finished_i = n.array(glob.glob( "/data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd"+str(nccd)+"_evt.fits" ))
	finished = n.array([ int(el.split('/')[-2]) for el in finished_i ])
	index.append( n.in1d(all_fields, finished) )

index = n.array((index))
sum_index = n.sum(index, axis=0)
completed_clu = (sum_index==7)

sum_index = n.sum([completed_clu, completed_agn, completed_bkg], axis=0)
completed_all = (sum_index==3)

already_processed_raw = n.array(glob.glob("/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/???/evt_???.fits"))
already_processed = n.array([ int(el.split('/')[-2]) for el in already_processed_raw])
to_keep = n.in1d(all_fields, already_processed, invert=True)

completed_not_processed = (completed_all) & (to_keep)

esass_code_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass')
HIDs, RAs, DECs = n.loadtxt(os.path.join(esass_code_dir, 'healpix_radius_nested_all.dat'), unpack=True)

f=open(os.path.join(esass_code_dir, 'calibrate_events_ready_to_go.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for HID, RA, DEC in zip(HIDs[completed_not_processed], RAs[completed_not_processed], DECs[completed_not_processed]):
	cmd = "python $GIT_AGN_MOCK/python/esass/calibrate_events.py "+str(int(HID)).zfill(3)+" "+str(RA)+" "+str(DEC)
	print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()	
