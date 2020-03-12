#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 613 
python ../pre-process-esass.py 613 
python ../esass.py 613 
