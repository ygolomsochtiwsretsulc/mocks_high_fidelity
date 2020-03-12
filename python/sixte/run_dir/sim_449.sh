#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 449 
python ../pre-process-esass.py 449 
python ../esass.py 449 
