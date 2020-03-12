#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 526 
python ../pre-process-esass.py 526 
python ../esass.py 526 
