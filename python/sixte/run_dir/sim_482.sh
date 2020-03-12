#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 482 
python ../pre-process-esass.py 482 
python ../esass.py 482 
