#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 745 
python ../pre-process-esass.py 745 
python ../esass.py 745 
