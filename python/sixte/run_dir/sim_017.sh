#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 017 
python ../pre-process-esass.py 017 
python ../esass.py 017 
