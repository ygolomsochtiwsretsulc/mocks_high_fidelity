import os, sys
import numpy as n

esass_code_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass')
esass_dir ="/data40s/erosim/eRASS/eSASS-v12-erass8-jc-clu-agn-bg/"
HIDs, RAs, DECs = n.loadtxt(os.path.join(esass_code_dir, 'healpix_radius_nested_all.dat'), unpack=True)

f=open(os.path.join(esass_code_dir, 'esass_commands_images.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for HID in n.arange(768):
	h_str = str(int(HID)).zfill(3)
	part1 = "nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_images.py "
	part2 = esass_dir+h_str+"/evt_"+h_str+".fits "
	part3 = esass_dir+h_str+"/ -b 16 --hp8"
	part4 = " > $GIT_AGN_MOCK/python/esass/logs/esass_image_b16_"+h_str+"_eSASSusers_200302.log & " 
	cmd = part1 + part2 + part3 + part4
	#print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()

f=open(os.path.join(esass_code_dir, 'esass_commands_detection.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for HID in n.arange(768):
	h_str = str(int(HID)).zfill(3)
	part1 = "nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py "
	part2 = esass_dir+h_str+" "
	part3 = esass_dir+h_str+" -V 13"
	part4 = " > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_"+h_str+"_eSASSusers_200302.log & " 
	cmd = part1 + part2 + part3 + part4
	#print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()

