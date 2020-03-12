#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 671 
python ../pre-process-esass.py 671 
python ../esass.py 671 
