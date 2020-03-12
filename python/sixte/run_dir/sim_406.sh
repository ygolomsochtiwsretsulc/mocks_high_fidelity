#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 406 
python ../pre-process-esass.py 406 
python ../esass.py 406 
