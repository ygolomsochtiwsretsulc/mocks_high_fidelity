#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 154 
python ../pre-process-esass.py 154 
python ../esass.py 154 
