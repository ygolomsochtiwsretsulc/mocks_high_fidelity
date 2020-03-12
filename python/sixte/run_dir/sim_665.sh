#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 665 
python ../pre-process-esass.py 665 
python ../esass.py 665 
