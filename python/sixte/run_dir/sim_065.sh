#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 065 
python ../pre-process-esass.py 065 
python ../esass.py 065 
