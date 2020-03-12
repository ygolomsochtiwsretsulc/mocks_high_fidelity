#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 464 
python ../pre-process-esass.py 464 
python ../esass.py 464 
