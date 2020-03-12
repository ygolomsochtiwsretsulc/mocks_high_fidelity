#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 467 
python ../pre-process-esass.py 467 
python ../esass.py 467 
