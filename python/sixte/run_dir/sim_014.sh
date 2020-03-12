#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 014 
python ../pre-process-esass.py 014 
python ../esass.py 014 
