#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 434 
python ../pre-process-esass.py 434 
python ../esass.py 434 
