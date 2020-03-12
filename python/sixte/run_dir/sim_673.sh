#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 673 
python ../pre-process-esass.py 673 
python ../esass.py 673 
