#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 043 
python ../pre-process-esass.py 043 
python ../esass.py 043 
