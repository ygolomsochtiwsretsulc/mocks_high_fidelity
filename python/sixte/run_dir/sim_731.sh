#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 731 
python ../pre-process-esass.py 731 
python ../esass.py 731 
