#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 747 
python ../pre-process-esass.py 747 
python ../esass.py 747 
