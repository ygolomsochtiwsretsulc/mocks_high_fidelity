#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 596 
python ../pre-process-esass.py 596 
python ../esass.py 596 
