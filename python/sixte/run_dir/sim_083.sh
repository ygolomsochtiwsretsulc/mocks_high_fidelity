#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 083 
python ../pre-process-esass.py 083 
python ../esass.py 083 
