#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 133 
python ../pre-process-esass.py 133 
python ../esass.py 133 
