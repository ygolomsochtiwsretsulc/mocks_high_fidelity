#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 495 
python ../pre-process-esass.py 495 
python ../esass.py 495 
