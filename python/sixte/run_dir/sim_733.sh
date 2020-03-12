#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 733 
python ../pre-process-esass.py 733 
python ../esass.py 733 
