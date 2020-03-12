#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 545 
python ../pre-process-esass.py 545 
python ../esass.py 545 
