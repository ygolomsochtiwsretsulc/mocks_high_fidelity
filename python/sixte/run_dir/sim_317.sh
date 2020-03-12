#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 317 
python ../pre-process-esass.py 317 
python ../esass.py 317 
