#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 724 
python ../pre-process-esass.py 724 
python ../esass.py 724 
