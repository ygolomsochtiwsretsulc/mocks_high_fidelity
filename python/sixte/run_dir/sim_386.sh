#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 386 
python ../pre-process-esass.py 386 
python ../esass.py 386 
