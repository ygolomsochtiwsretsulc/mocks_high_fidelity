#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 391 
python ../pre-process-esass.py 391 
python ../esass.py 391 
