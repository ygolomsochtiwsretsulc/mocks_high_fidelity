#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 387 
python ../pre-process-esass.py 387 
python ../esass.py 387 
