#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 330 
python ../pre-process-esass.py 330 
python ../esass.py 330 
