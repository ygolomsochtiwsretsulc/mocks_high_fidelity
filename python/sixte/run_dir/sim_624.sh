#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 624 
python ../pre-process-esass.py 624 
python ../esass.py 624 
