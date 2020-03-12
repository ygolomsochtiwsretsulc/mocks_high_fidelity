#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 011 
python ../pre-process-esass.py 011 
python ../esass.py 011 
