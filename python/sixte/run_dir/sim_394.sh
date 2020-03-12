#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 394 
python ../pre-process-esass.py 394 
python ../esass.py 394 
