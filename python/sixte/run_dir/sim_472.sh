#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 472 
python ../pre-process-esass.py 472 
python ../esass.py 472 
