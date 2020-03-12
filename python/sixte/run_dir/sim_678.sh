#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 678 
python ../pre-process-esass.py 678 
python ../esass.py 678 
