#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 759 
python ../pre-process-esass.py 759 
python ../esass.py 759 
