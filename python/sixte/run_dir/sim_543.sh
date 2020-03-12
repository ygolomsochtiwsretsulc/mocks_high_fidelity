#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 543 
python ../pre-process-esass.py 543 
python ../esass.py 543 
