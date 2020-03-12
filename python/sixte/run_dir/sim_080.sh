#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 080 
python ../pre-process-esass.py 080 
python ../esass.py 080 
