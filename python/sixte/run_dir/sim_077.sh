#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 077 
python ../pre-process-esass.py 077 
python ../esass.py 077 
