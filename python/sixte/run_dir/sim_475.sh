#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 475 
python ../pre-process-esass.py 475 
python ../esass.py 475 
