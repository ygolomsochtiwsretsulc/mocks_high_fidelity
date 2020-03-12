#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 337 
python ../pre-process-esass.py 337 
python ../esass.py 337 
