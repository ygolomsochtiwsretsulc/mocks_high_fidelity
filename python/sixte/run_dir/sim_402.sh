#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 402 
python ../pre-process-esass.py 402 
python ../esass.py 402 
