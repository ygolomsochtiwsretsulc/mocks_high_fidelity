#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 766 
python ../pre-process-esass.py 766 
python ../esass.py 766 
