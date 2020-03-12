#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 555 
python ../pre-process-esass.py 555 
python ../esass.py 555 
