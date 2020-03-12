#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 636 
python ../pre-process-esass.py 636 
python ../esass.py 636 
