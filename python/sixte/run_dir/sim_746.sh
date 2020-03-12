#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 746 
python ../pre-process-esass.py 746 
python ../esass.py 746 
