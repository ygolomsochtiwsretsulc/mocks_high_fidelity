#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 705 
python ../pre-process-esass.py 705 
python ../esass.py 705 
