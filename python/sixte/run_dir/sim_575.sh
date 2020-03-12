#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 575 
python ../pre-process-esass.py 575 
python ../esass.py 575 
