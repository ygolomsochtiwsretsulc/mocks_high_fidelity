#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 557 
python ../pre-process-esass.py 557 
python ../esass.py 557 
