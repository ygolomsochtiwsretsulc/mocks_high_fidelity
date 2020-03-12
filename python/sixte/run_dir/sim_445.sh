#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 445 
python ../pre-process-esass.py 445 
python ../esass.py 445 
