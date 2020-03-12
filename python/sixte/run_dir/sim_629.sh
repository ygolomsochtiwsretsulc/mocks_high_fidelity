#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 629 
python ../pre-process-esass.py 629 
python ../esass.py 629 
