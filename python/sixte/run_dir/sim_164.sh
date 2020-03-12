#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 164 
python ../pre-process-esass.py 164 
python ../esass.py 164 
