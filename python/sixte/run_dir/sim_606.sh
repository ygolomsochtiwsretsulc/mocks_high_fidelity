#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 606 
python ../pre-process-esass.py 606 
python ../esass.py 606 
