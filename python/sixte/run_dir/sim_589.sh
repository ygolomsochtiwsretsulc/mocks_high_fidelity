#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 589 
python ../pre-process-esass.py 589 
python ../esass.py 589 
