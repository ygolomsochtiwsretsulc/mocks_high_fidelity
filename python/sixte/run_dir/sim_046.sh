#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 046 
python ../pre-process-esass.py 046 
python ../esass.py 046 
