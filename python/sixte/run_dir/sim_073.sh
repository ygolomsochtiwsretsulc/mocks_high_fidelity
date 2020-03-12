#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 073 
python ../pre-process-esass.py 073 
python ../esass.py 073 
