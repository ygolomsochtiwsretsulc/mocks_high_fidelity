#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 651 
python ../pre-process-esass.py 651 
python ../esass.py 651 
