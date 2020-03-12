#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 201 
python ../pre-process-esass.py 201 
python ../esass.py 201 
