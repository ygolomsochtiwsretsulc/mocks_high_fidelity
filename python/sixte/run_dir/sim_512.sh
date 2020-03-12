#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 512 
python ../pre-process-esass.py 512 
python ../esass.py 512 
