#!/bin/bash 
source /home/erosita/sw/sass-setup.sh eSASSdevel 
python ../simulate.py 090 
python ../pre-process-esass.py 090 
python ../esass.py 090 
