#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 676 
python ../pre-process-esass.py 676 
python ../esass.py 676 
