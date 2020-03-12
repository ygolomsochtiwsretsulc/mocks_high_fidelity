#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 428 
python ../pre-process-esass.py 428 
python ../esass.py 428 
