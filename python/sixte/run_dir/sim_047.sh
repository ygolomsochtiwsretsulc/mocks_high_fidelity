#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 047 
python ../pre-process-esass.py 047 
python ../esass.py 047 
