#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 036 
python ../pre-process-esass.py 036 
python ../esass.py 036 
