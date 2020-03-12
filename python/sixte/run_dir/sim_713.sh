#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 713 
python ../pre-process-esass.py 713 
python ../esass.py 713 
