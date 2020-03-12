#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 566 
python ../pre-process-esass.py 566 
python ../esass.py 566 
