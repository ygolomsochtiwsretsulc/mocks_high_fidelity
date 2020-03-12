#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 698 
python ../pre-process-esass.py 698 
python ../esass.py 698 
