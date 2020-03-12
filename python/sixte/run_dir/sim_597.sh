#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 597 
python ../pre-process-esass.py 597 
python ../esass.py 597 
