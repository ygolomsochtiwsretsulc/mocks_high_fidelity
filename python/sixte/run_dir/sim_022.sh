#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 022 
python ../pre-process-esass.py 022 
python ../esass.py 022 
