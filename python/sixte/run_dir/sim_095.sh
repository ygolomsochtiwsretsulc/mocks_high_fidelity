#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 095 
python ../pre-process-esass.py 095 
python ../esass.py 095 
