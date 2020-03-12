#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 418 
python ../pre-process-esass.py 418 
python ../esass.py 418 
