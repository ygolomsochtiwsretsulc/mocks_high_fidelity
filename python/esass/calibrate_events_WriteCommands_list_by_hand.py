import os, sys
import numpy as n

esass_code_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass')
HIDs, RAs, DECs = n.loadtxt(os.path.join(esass_code_dir, 'healpix_radius_nested_all.dat'), unpack=True)

HID_selection = n.array([ 157, 
 180, 
 181, 
 182, 
 183, 
 238, 
 250, 
 354, 
 383, 
 387, 
 517, 
 586, 
 675, 
 684, 
 689, 
 691, 
 704, 
 713, 
 721, 
 724, 
 735
 ])

f=open(os.path.join(esass_code_dir, 'calibrate_events_list_by_hand.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for HID, RA, DEC in zip(HIDs[HID_selection], RAs[HID_selection], DECs[HID_selection]):
	cmd = "python $GIT_AGN_MOCK/python/esass/calibrate_events.py "+str(int(HID)).zfill(3)+" "+str(RA)+" "+str(DEC)
	print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()	
