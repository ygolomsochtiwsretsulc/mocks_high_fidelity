#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 567 
python ../pre-process-esass.py 567 
python ../esass.py 567 
