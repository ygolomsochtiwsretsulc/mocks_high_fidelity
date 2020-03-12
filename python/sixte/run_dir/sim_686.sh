#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 686 
python ../pre-process-esass.py 686 
python ../esass.py 686 
