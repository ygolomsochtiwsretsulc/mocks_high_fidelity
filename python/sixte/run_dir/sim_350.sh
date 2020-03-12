#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 350 
python ../pre-process-esass.py 350 
python ../esass.py 350 
