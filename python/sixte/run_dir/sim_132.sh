#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 132 
python ../pre-process-esass.py 132 
python ../esass.py 132 
