#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 019 
python ../pre-process-esass.py 019 
python ../esass.py 019 
