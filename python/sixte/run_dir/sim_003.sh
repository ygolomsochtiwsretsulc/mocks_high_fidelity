#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 003 
python ../pre-process-esass.py 003 
python ../esass.py 003 
