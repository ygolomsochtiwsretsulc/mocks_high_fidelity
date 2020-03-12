#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 757 
python ../pre-process-esass.py 757 
python ../esass.py 757 
