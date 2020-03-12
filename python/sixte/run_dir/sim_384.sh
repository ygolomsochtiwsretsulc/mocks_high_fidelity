#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 384 
python ../pre-process-esass.py 384 
python ../esass.py 384 
