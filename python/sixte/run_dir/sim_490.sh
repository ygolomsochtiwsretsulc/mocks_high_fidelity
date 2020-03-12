#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 490 
python ../pre-process-esass.py 490 
python ../esass.py 490 
