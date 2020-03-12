#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 681 
python ../pre-process-esass.py 681 
python ../esass.py 681 
