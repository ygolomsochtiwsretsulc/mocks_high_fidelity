#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 579 
python ../pre-process-esass.py 579 
python ../esass.py 579 
