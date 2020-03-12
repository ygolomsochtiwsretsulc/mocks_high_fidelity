#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 646 
python ../pre-process-esass.py 646 
python ../esass.py 646 
