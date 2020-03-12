#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 435 
python ../pre-process-esass.py 435 
python ../esass.py 435 
