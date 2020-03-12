#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 484 
python ../pre-process-esass.py 484 
python ../esass.py 484 
