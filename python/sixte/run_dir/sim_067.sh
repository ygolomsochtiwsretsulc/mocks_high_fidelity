#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 067 
python ../pre-process-esass.py 067 
python ../esass.py 067 
