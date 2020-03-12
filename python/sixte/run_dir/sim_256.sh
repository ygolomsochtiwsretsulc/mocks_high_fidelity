#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 256 
python ../pre-process-esass.py 256 
python ../esass.py 256 
