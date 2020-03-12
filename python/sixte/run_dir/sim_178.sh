#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 178 
python ../pre-process-esass.py 178 
python ../esass.py 178 
