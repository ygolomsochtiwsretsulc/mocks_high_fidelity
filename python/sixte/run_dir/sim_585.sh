#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 585 
python ../pre-process-esass.py 585 
python ../esass.py 585 
