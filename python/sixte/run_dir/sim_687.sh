#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 687 
python ../pre-process-esass.py 687 
python ../esass.py 687 
