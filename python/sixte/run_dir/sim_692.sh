#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 692 
python ../pre-process-esass.py 692 
python ../esass.py 692 
