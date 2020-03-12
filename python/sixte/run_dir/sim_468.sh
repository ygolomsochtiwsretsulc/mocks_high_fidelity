#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 468 
python ../pre-process-esass.py 468 
python ../esass.py 468 
