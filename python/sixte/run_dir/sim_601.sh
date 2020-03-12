#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 601 
python ../pre-process-esass.py 601 
python ../esass.py 601 
