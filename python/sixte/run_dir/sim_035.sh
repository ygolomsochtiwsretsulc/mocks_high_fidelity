#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 035 
python ../pre-process-esass.py 035 
python ../esass.py 035 
