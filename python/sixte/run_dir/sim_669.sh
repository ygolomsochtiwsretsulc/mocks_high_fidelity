#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 669 
python ../pre-process-esass.py 669 
python ../esass.py 669 
