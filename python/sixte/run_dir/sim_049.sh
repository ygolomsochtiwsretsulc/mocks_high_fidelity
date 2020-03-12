#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 049 
python ../pre-process-esass.py 049 
python ../esass.py 049 
