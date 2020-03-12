#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 548 
python ../pre-process-esass.py 548 
python ../esass.py 548 
