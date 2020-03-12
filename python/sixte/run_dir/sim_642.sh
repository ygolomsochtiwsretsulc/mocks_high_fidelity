#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 642 
python ../pre-process-esass.py 642 
python ../esass.py 642 
