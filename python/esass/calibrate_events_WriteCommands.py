import os, sys
import numpy as n

esass_code_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass')
HIDs, RAs, DECs = n.loadtxt(os.path.join(esass_code_dir, 'healpix_radius_nested_all.dat'), unpack=True)

f=open(os.path.join(esass_code_dir, 'calibrate_events.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for HID, RA, DEC in zip(HIDs, RAs, DECs):
	cmd = "python $GIT_AGN_MOCK/python/esass/calibrate_events.py "+str(int(HID)).zfill(3)+" "+str(RA)+" "+str(DEC)
	print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()	
