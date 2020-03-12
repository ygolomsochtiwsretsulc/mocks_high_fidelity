#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 358 
python ../pre-process-esass.py 358 
python ../esass.py 358 
