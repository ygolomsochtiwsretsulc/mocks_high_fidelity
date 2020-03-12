#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 699 
python ../pre-process-esass.py 699 
python ../esass.py 699 
