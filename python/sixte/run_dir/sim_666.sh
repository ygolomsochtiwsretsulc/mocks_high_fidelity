#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 666 
python ../pre-process-esass.py 666 
python ../esass.py 666 
