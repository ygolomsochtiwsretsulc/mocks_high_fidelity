#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 117 
python ../pre-process-esass.py 117 
python ../esass.py 117 
