#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 060 
python ../pre-process-esass.py 060 
python ../esass.py 060 
