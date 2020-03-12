#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 718 
python ../pre-process-esass.py 718 
python ../esass.py 718 
