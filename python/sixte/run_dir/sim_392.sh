#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 392 
python ../pre-process-esass.py 392 
python ../esass.py 392 
