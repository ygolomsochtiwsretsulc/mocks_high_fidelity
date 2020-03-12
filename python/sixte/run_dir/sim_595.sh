#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 595 
python ../pre-process-esass.py 595 
python ../esass.py 595 
