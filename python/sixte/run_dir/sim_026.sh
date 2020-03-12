#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 026 
python ../pre-process-esass.py 026 
python ../esass.py 026 
