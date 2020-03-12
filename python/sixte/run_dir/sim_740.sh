#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 740 
python ../pre-process-esass.py 740 
python ../esass.py 740 
