#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 755 
python ../pre-process-esass.py 755 
python ../esass.py 755 
