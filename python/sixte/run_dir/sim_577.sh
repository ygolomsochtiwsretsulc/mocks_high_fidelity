#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 577 
python ../pre-process-esass.py 577 
python ../esass.py 577 
