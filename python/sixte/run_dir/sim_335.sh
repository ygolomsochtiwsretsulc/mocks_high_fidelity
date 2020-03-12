#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 335 
python ../pre-process-esass.py 335 
python ../esass.py 335 
