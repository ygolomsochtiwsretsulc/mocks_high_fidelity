#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 282 
python ../pre-process-esass.py 282 
python ../esass.py 282 
