#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 568 
python ../pre-process-esass.py 568 
python ../esass.py 568 
