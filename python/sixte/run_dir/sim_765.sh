#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 765 
python ../pre-process-esass.py 765 
python ../esass.py 765 
