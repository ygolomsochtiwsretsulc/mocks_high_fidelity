#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 424 
python ../pre-process-esass.py 424 
python ../esass.py 424 
