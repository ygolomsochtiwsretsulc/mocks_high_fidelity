#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 355 
python ../pre-process-esass.py 355 
python ../esass.py 355 
