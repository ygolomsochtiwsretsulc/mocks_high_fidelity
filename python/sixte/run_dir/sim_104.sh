#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 104 
python ../pre-process-esass.py 104 
python ../esass.py 104 
