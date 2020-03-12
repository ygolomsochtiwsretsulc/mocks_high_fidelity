#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 576 
python ../pre-process-esass.py 576 
python ../esass.py 576 
