#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 368 
python ../pre-process-esass.py 368 
python ../esass.py 368 
