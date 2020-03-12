#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 610 
python ../pre-process-esass.py 610 
python ../esass.py 610 
