#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 533 
python ../pre-process-esass.py 533 
python ../esass.py 533 
