#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 225 
python ../pre-process-esass.py 225 
python ../esass.py 225 
