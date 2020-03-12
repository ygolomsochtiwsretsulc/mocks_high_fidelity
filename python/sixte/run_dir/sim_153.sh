#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 153 
python ../pre-process-esass.py 153 
python ../esass.py 153 
