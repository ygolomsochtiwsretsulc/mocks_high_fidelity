#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 647 
python ../pre-process-esass.py 647 
python ../esass.py 647 
