#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 761 
python ../pre-process-esass.py 761 
python ../esass.py 761 
