#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 184 
python ../pre-process-esass.py 184 
python ../esass.py 184 
