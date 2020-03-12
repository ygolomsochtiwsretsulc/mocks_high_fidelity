#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 031 
python ../pre-process-esass.py 031 
python ../esass.py 031 
