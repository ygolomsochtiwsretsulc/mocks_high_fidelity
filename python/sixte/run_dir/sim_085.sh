#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 085 
python ../pre-process-esass.py 085 
python ../esass.py 085 
