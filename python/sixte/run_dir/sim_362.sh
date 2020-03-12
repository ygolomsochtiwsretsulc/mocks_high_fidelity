#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 362 
python ../pre-process-esass.py 362 
python ../esass.py 362 
