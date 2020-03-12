#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 082 
python ../pre-process-esass.py 082 
python ../esass.py 082 
