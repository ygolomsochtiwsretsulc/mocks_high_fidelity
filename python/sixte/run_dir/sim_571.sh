#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 571 
python ../pre-process-esass.py 571 
python ../esass.py 571 
