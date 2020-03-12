#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 726 
python ../pre-process-esass.py 726 
python ../esass.py 726 
