#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 025 
python ../pre-process-esass.py 025 
python ../esass.py 025 
