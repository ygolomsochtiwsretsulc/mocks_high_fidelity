#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 272 
python ../pre-process-esass.py 272 
python ../esass.py 272 
