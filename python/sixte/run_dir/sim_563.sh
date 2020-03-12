#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 563 
python ../pre-process-esass.py 563 
python ../esass.py 563 
