#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 691 
python ../pre-process-esass.py 691 
python ../esass.py 691 
