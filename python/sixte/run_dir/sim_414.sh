#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 414 
python ../pre-process-esass.py 414 
python ../esass.py 414 
