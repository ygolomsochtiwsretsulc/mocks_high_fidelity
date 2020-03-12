#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 397 
python ../pre-process-esass.py 397 
python ../esass.py 397 
