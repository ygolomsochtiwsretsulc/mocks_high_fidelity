#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 721 
python ../pre-process-esass.py 721 
python ../esass.py 721 
