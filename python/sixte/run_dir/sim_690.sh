#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 690 
python ../pre-process-esass.py 690 
python ../esass.py 690 
