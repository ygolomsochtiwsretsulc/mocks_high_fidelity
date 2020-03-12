#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 668 
python ../pre-process-esass.py 668 
python ../esass.py 668 
