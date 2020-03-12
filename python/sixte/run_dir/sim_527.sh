#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 527 
python ../pre-process-esass.py 527 
python ../esass.py 527 
