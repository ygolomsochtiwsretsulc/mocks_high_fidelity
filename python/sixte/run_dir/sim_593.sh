#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 593 
python ../pre-process-esass.py 593 
python ../esass.py 593 
