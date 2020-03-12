#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 315 
python ../pre-process-esass.py 315 
python ../esass.py 315 
