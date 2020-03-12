#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 588 
python ../pre-process-esass.py 588 
python ../esass.py 588 
