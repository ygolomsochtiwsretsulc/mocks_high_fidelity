#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 257 
python ../pre-process-esass.py 257 
python ../esass.py 257 
