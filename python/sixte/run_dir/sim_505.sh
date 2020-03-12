#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 505 
python ../pre-process-esass.py 505 
python ../esass.py 505 
