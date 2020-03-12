#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 515 
python ../pre-process-esass.py 515 
python ../esass.py 515 
