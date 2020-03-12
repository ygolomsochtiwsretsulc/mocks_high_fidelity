#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 663 
python ../pre-process-esass.py 663 
python ../esass.py 663 
