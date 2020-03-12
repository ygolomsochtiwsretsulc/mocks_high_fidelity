#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 103 
python ../pre-process-esass.py 103 
python ../esass.py 103 
