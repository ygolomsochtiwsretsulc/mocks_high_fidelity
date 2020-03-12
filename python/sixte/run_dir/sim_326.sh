#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 326 
python ../pre-process-esass.py 326 
python ../esass.py 326 
