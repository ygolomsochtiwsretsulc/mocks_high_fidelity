#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 000 
python ../pre-process-esass.py 000 
python ../esass.py 000 
