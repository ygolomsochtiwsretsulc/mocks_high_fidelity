#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 583 
python ../pre-process-esass.py 583 
python ../esass.py 583 
