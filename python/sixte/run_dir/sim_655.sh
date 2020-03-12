#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 655 
python ../pre-process-esass.py 655 
python ../esass.py 655 
