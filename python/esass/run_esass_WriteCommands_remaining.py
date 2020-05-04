import numpy as n
import glob, os, sys
code_dir = os.path.join(os.environ['GIT_AGN_MOCK'], 'python', 'esass')
esass_dir ="/data26s/comparat/simulations/erosim/eSASS-v12-erass8-jc-clu-agn-bg/"
all_fields = n.arange(768)
# determine all_fields not in all_fields
index = []
for nccd in n.arange(7)+1:
	finished_i = n.array(glob.glob( "/data26s/comparat/simulations/erosim/eSASS-v12-erass8-jc-clu-agn-bg/???/???_0216_DetMsk.fits" ))
	finished = n.array([ int(el.split('/')[-2]) for el in finished_i ])
	index.append( n.in1d(all_fields, finished) )

index = n.array((index))
sum_index = n.sum(index, axis=0)
completed = (sum_index==7)
remaining = (completed==False)
completed_image = n.copy(completed)

f=open(os.path.join(code_dir, 'esass_commands_images_remaining.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for el in all_fields[remaining]:
	h_str = str(el).zfill(3)
	part1 = "nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_images.py "
	part2 = esass_dir+h_str+"/evt_"+h_str+".fits "
	part3 = esass_dir+h_str+"/ -b 16 --hp8"
	part4 = " > $GIT_AGN_MOCK/python/esass/logs/esass_image_b16_"+h_str+"_eSASSusers_200302.log & " 
	cmd = part1 + part2 + part3 + part4
	print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()

all_fields = n.arange(768)
# determine all_fields not in all_fields
index = []
for nccd in n.arange(7)+1:
	finished_i = n.array(glob.glob( "/data26s/comparat/simulations/erosim/eSASS-v12-erass8-jc-clu-agn-bg/???/???_0216_Sc1Cat.fits" ))
	finished = n.array([ int(el.split('/')[-2]) for el in finished_i ])
	index.append( n.in1d(all_fields, finished) )

index = n.array((index))
sum_index = n.sum(index, axis=0)
completed = (sum_index==7)
remaining = (completed==False)&(completed_image)

f=open(os.path.join(code_dir, 'esass_commands_detection_remaining.sh'), 'w')
f.write("""#!/bin/bash/""")
f.write('\n')
for el in all_fields[remaining]:
	h_str = str(el).zfill(3)
	part1 = "nohup python $GIT_AGN_MOCK/python/esass/eRO_pipe_det.py "
	part2 = esass_dir+h_str+" "
	part3 = esass_dir+h_str+" -V 13"
	part4 = " > $GIT_AGN_MOCK/python/esass/logs/esass_detpipe_v13_"+h_str+"_eSASSusers_200302.log & " 
	cmd = part1 + part2 + part3 + part4
	print(cmd)
	f.write(cmd)
	f.write('\n')
	
f.close()

