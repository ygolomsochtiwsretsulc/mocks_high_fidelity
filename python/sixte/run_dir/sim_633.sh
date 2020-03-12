#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 633 
python ../pre-process-esass.py 633 
python ../esass.py 633 
