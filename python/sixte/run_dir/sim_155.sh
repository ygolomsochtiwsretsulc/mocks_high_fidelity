#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 155 
python ../pre-process-esass.py 155 
python ../esass.py 155 
