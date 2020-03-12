#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 214 
python ../pre-process-esass.py 214 
python ../esass.py 214 
