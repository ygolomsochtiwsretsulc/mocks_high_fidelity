#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 229 
python ../pre-process-esass.py 229 
python ../esass.py 229 
