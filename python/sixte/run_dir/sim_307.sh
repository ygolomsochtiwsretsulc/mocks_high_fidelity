#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 307 
python ../pre-process-esass.py 307 
python ../esass.py 307 
