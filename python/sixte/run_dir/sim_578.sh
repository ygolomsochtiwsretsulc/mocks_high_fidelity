#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 578 
python ../pre-process-esass.py 578 
python ../esass.py 578 
