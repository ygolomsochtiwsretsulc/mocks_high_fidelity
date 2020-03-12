#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 021 
python ../pre-process-esass.py 021 
python ../esass.py 021 
