#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 564 
python ../pre-process-esass.py 564 
python ../esass.py 564 
