#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 278 
python ../pre-process-esass.py 278 
python ../esass.py 278 
