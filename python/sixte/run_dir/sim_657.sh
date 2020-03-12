#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 657 
python ../pre-process-esass.py 657 
python ../esass.py 657 
