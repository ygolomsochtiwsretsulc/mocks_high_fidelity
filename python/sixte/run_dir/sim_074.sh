#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 074 
python ../pre-process-esass.py 074 
python ../esass.py 074 
