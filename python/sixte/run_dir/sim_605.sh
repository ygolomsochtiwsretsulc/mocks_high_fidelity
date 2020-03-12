#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 605 
python ../pre-process-esass.py 605 
python ../esass.py 605 
