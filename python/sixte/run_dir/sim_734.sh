#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 734 
python ../pre-process-esass.py 734 
python ../esass.py 734 
