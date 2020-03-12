#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 440 
python ../pre-process-esass.py 440 
python ../esass.py 440 
