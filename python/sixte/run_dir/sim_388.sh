#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 388 
python ../pre-process-esass.py 388 
python ../esass.py 388 
