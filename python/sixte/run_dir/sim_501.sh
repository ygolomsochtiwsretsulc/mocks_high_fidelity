#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 501 
python ../pre-process-esass.py 501 
python ../esass.py 501 
