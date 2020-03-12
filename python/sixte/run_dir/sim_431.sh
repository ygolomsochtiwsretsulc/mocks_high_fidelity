#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 431 
python ../pre-process-esass.py 431 
python ../esass.py 431 
