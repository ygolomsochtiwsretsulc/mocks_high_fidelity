#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 607 
python ../pre-process-esass.py 607 
python ../esass.py 607 
