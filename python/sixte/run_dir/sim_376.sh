#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 376 
python ../pre-process-esass.py 376 
python ../esass.py 376 
