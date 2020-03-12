#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 478 
python ../pre-process-esass.py 478 
python ../esass.py 478 
