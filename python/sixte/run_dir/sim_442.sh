#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 442 
python ../pre-process-esass.py 442 
python ../esass.py 442 
