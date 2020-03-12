#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 268 
python ../pre-process-esass.py 268 
python ../esass.py 268 
