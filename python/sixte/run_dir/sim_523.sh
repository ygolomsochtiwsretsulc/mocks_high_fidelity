#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 523 
python ../pre-process-esass.py 523 
python ../esass.py 523 
