#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 627 
python ../pre-process-esass.py 627 
python ../esass.py 627 
