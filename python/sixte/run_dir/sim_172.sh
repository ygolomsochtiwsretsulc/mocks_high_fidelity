#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 172 
python ../pre-process-esass.py 172 
python ../esass.py 172 
