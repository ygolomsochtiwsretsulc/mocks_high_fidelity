#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 626 
python ../pre-process-esass.py 626 
python ../esass.py 626 
