#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 451 
python ../pre-process-esass.py 451 
python ../esass.py 451 
