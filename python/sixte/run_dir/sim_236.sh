#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 236 
python ../pre-process-esass.py 236 
python ../esass.py 236 
