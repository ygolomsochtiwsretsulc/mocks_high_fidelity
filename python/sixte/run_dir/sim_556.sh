#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 556 
python ../pre-process-esass.py 556 
python ../esass.py 556 
