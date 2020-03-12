#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 659 
python ../pre-process-esass.py 659 
python ../esass.py 659 
