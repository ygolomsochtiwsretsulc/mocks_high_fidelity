#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 739 
python ../pre-process-esass.py 739 
python ../esass.py 739 
