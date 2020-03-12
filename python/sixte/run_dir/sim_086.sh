#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 086 
python ../pre-process-esass.py 086 
python ../esass.py 086 
