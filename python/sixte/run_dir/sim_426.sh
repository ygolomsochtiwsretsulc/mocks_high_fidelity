#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 426 
python ../pre-process-esass.py 426 
python ../esass.py 426 
