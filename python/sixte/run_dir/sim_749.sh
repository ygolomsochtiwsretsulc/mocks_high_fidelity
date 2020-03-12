#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 749 
python ../pre-process-esass.py 749 
python ../esass.py 749 
