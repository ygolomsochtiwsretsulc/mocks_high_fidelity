#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 007 
python ../pre-process-esass.py 007 
python ../esass.py 007 
