#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 722 
python ../pre-process-esass.py 722 
python ../esass.py 722 
