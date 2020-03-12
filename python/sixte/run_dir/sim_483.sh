#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 483 
python ../pre-process-esass.py 483 
python ../esass.py 483 
