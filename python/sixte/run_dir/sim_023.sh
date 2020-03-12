#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 023 
python ../pre-process-esass.py 023 
python ../esass.py 023 
