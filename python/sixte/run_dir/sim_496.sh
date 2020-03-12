#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 496 
python ../pre-process-esass.py 496 
python ../esass.py 496 
