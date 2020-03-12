#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 033 
python ../pre-process-esass.py 033 
python ../esass.py 033 
