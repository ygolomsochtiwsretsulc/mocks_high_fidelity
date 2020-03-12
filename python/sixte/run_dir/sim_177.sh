#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 177 
python ../pre-process-esass.py 177 
python ../esass.py 177 
