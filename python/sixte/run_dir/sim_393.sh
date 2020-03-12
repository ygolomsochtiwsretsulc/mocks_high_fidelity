#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 393 
python ../pre-process-esass.py 393 
python ../esass.py 393 
