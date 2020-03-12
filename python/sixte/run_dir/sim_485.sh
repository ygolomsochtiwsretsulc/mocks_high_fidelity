#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 485 
python ../pre-process-esass.py 485 
python ../esass.py 485 
