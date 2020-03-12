#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 704 
python ../pre-process-esass.py 704 
python ../esass.py 704 
