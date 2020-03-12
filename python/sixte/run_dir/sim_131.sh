#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 131 
python ../pre-process-esass.py 131 
python ../esass.py 131 
