#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 068 
python ../pre-process-esass.py 068 
python ../esass.py 068 
