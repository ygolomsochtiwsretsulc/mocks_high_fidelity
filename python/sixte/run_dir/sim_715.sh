#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 715 
python ../pre-process-esass.py 715 
python ../esass.py 715 
