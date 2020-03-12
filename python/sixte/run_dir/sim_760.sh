#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 760 
python ../pre-process-esass.py 760 
python ../esass.py 760 
