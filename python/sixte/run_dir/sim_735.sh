#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 735 
python ../pre-process-esass.py 735 
python ../esass.py 735 
