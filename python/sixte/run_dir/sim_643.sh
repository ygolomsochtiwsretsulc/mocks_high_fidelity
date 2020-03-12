#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 643 
python ../pre-process-esass.py 643 
python ../esass.py 643 
