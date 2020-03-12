#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 685 
python ../pre-process-esass.py 685 
python ../esass.py 685 
