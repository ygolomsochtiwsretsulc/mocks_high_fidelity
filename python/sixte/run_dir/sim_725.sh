#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 725 
python ../pre-process-esass.py 725 
python ../esass.py 725 
