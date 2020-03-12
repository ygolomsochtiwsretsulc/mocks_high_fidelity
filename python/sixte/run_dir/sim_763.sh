#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 763 
python ../pre-process-esass.py 763 
python ../esass.py 763 
