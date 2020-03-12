#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 614 
python ../pre-process-esass.py 614 
python ../esass.py 614 
