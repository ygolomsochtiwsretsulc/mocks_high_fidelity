#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 616 
python ../pre-process-esass.py 616 
python ../esass.py 616 
