#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 267 
python ../pre-process-esass.py 267 
python ../esass.py 267 
