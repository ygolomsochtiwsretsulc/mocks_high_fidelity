#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 158 
python ../pre-process-esass.py 158 
python ../esass.py 158 
