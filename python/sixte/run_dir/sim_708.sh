#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 708 
python ../pre-process-esass.py 708 
python ../esass.py 708 
