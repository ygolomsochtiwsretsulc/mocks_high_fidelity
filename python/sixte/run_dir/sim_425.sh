#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 425 
python ../pre-process-esass.py 425 
python ../esass.py 425 
