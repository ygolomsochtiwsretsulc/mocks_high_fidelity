#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 473 
python ../pre-process-esass.py 473 
python ../esass.py 473 
