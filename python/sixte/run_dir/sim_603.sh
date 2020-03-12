#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 603 
python ../pre-process-esass.py 603 
python ../esass.py 603 
