#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 075 
python ../pre-process-esass.py 075 
python ../esass.py 075 
