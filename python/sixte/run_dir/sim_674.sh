#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 674 
python ../pre-process-esass.py 674 
python ../esass.py 674 
