#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 302 
python ../pre-process-esass.py 302 
python ../esass.py 302 
