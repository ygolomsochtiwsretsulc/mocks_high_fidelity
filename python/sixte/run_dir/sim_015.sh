#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 015 
python ../pre-process-esass.py 015 
python ../esass.py 015 
