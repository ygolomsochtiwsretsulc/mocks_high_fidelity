#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 040 
python ../pre-process-esass.py 040 
python ../esass.py 040 
