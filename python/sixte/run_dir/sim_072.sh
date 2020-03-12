#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 072 
python ../pre-process-esass.py 072 
python ../esass.py 072 
