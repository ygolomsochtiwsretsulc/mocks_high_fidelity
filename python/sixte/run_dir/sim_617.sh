#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 617 
python ../pre-process-esass.py 617 
python ../esass.py 617 
