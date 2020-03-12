#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 071 
python ../pre-process-esass.py 071 
python ../esass.py 071 
